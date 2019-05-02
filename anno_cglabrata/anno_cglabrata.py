from sequence_lib import rc_seq
from sequence_lib import translate_exon
from biofile import simple_fasta_write
from biofile import simple_fasta_load
import os
import re
import sys
import copy
import re


def annotate(args):
    orf_calling(
        fa = args['query'],
        cdna= args['ref'],
        align=args['align'],
        min_length=args['min_length'],
        prefix=args['prefix'],
        strain=args['strain']
    )


def orf_calling(fa: str, cdna: str, align: str, min_length: int, prefix: str, strain='BG2'):
    chroms, seqs = simple_fasta_load(fa)
    seqs = [s.upper() for s in seqs]
    genome = zip(chroms, seqs)
    refgenes, refseqs = simple_fasta_load(cdna)
    refgenes.sort()

    ecod = ['TAG', 'TGA', 'TAA']
    scod = 'ATG'
    rscod = 'CAT'
    recod = ['CTA', 'TCA', 'TTA']
    orf_loc = {}

    # Load alignment
    best_orf_align = {}
    orf_align = {}
    for c in chroms:
        orf_align[c] = {i: [] for i in range(6)}
    with open(align) as filep:
        for line in filep:
            ent = line.split()
            gene = ent[0]
            chrom = ent[1]
            gs = int(ent[6])
            ss = int(ent[8])
            se = int(ent[9])
            score = float(ent[-1])
            if int(ent[3]) < min_length:
                continue
            direction = 1 if ss < se else -1
            if direction == 1:
                frame = (ss - gs) % 3
                orf_align[chrom][frame].append([ss - 1, gene, score])
                if gene not in best_orf_align:
                    best_orf_align[gene] = [chrom, ss - 1]
            else:
                frame = (ss + gs - 1) % 3
                orf_align[chrom][frame + 3].append([ss, gene, score])
                if gene not in best_orf_align:
                    best_orf_align[gene] = [chrom, ss]
    print(len(best_orf_align))

    # Raw orf map
    for c, s in genome:
        # +3 +1 +2, -1, -2, -3
        orf_loc[c] = [{} for i in range(6)]
        # W strand
        start = [[m.start(), m.start() % 3] for m in re.finditer(scod, s)]
        end = []
        for e in ecod:
            end.extend([m.end(), m.end() % 3] for m in re.finditer(e, s))
        end.sort()
        # C strand
        rstart = [[m.end(), m.end() % 3] for m in re.finditer(rscod, s)]
        rstart.sort(reverse=True)
        rend = []
        for e in recod:
            rend.extend([m.start(), m.start() % 3] for m in re.finditer(e, s))
        rend.sort(reverse=True)
        for i in range(3):
            temp_start = [s[0] for s in start if s[1] == i]
            temp_end = [e[0] for e in end if e[1] == i]
            p_start = 0
            for j in range(len(temp_end)):
                cur_end = temp_end[j]
                cur_start = temp_start[p_start]
                if cur_end - cur_start > min_length:
                    orf_loc[c][i][(cur_start, cur_end)] = []
                while p_start < len(temp_start) - 1 and temp_start[p_start] < cur_end:
                    p_start += 1
        for i in range(3):
            temp_start = [s[0] for s in rstart if s[1] == i]
            temp_end = [e[0] for e in rend if e[1] == i]
            p_start = 0
            for j in range(len(temp_end)):
                cur_end = temp_end[j]
                cur_start = temp_start[p_start]
                if cur_start - cur_end > min_length:
                    orf_loc[c][3 + i][(cur_start, cur_end)] = []
                while p_start < len(temp_start) - 1 and temp_start[p_start] > cur_end:
                    p_start += 1

    # Load evidence
    for c, oloc in orf_loc.items():
        count = 0
        for f, loc in enumerate(oloc):
            aloc = orf_align[c][f]
            print(c, len(orf_loc[c][f]), 'Orf Prediction', len(aloc))
            for aln in aloc:
                for reg in loc.keys():
                    if f < 3:
                        if reg[0] <= aln[0] < reg[1]:
                            orf_loc[c][f][reg].append((aln[1], aln[2]))
                            break
                    else:
                        if reg[0] >= aln[0] > reg[1]:
                            orf_loc[c][f][reg].append((aln[1], aln[2]))
                            break

    # supported orfs
    sup_orf_loc = {}
    for c, oloc in orf_loc.items():
        sup_orf_loc[c] = {}
        count = 0
        for f, loc in enumerate(oloc):
            for l, evi in loc.items():
                if evi:
                    evi.sort(key=lambda x: (x[0], -x[1]))
                    combine_evi = []
                    cur_ev = list(evi[0])
                    for e in evi[1:]:
                        if e[0] == cur_ev[0]:
                            cur_ev[1] += e[1]
                            # continue
                        else:
                            combine_evi.append(tuple(cur_ev))
                            cur_ev = list(e)
                    combine_evi.append(tuple(cur_ev))
                    if f < 3:
                        sup_orf_loc[c][(l[0], l[1], 1)] = combine_evi
                    else:
                        sup_orf_loc[c][(l[1], l[0], -1)] = combine_evi
                    count += len(combine_evi)
        print('Supported Evi: %s: %d: %d' % (c, len(sup_orf_loc[c]), count))

    non_overlap = {}
    non_overlap_id = {}

    non_overlap = {}
    # Remove totally overalped orfs
    for c in sup_orf_loc:
        non_overlap[c] = {}
        locs = sorted(list(sup_orf_loc[c].keys()), key=lambda x: (x[0], -x[1]))
        if len(locs) == 0:
            continue
        else:
            cur_loc = locs[0]
            for l in locs[1:]:
                # Totally overlapped
                if cur_loc[0] <= l[0] <= l[1] <= cur_loc[1]:
                    continue
                else:
                    non_overlap[c][cur_loc] = sup_orf_loc[c][cur_loc]
                    cur_loc = l
            non_overlap[c][cur_loc] = sup_orf_loc[c][cur_loc]

    total_non_overlap = 0
    for c in non_overlap:
        regs = list(non_overlap[c].keys())
        regs.sort()
        total_non_overlap += len(regs)
        non_overlap_id[c] = {reg: i for i, reg in enumerate(regs)}

    # Reciprocal_best
    rb_assign = {}
    gene_rb = []
    gene_rb_assign = {}
    non_rb_assign = {}
    for c, loc in non_overlap.items():
        rb_assign[c] = {}
        non_rb_assign[c] = {}
        for reg, evi in loc.items():
            evi.sort(key=lambda x: x[1], reverse=True)
            best_gene = evi[0][0]
            gene_best_loc = best_orf_align[best_gene]
            if gene_best_loc[0] != c:
                non_rb_assign[c][reg] = evi
                continue
            if reg[0] <= gene_best_loc[1] <= reg[1] or reg[0] >= gene_best_loc[1] >= reg[1]:
                rb_assign[c][reg] = best_gene
                gene_rb.append(best_gene)
                gene_rb_assign[best_gene] = [c, reg[0], reg[1], reg[2]]
                continue
            else:
                non_rb_assign[c][reg] = evi

    gene_rb.sort()

    # Save recipocal best
    with open('%s.gene.reciprocal.tab' % prefix, 'w') as filep:
        filep.write('\n'.join(gene_rb))

    non_rb_loc = {}
    # Remove Rb gene in non-rb orfs
    for c, loc in non_rb_assign.items():
        for reg, evi in loc.items():
            new_evi = [e for e in evi if e[0] not in gene_rb]
            non_rb_assign[c][reg] = new_evi
            for e in new_evi:
                non_rb_loc.setdefault(e[0], [])
                non_rb_loc[e[0]].append((c, reg[0], reg[1], reg[2]))

    gene_syn_assign = {}
    gene_syn = {}
    # Easy synteny assign
    inter_chrom_penalty = len(refgenes) // 13
    non_rb_genes = list(non_rb_loc.keys())
    non_rb_genes.sort()
    for g in non_rb_genes:
        loc = non_rb_loc[g]
        g1, g2, _ = search_synteny(g, gene_rb)
        g1rid = refgenes.index(g1)
        g2rid = refgenes.index(g2)
        g1loc = gene_rb_assign[g1]
        g2loc = gene_rb_assign[g2]
        g1qid = non_overlap_id[g1loc[0]][tuple(g1loc[1:])]
        g2qid = non_overlap_id[g2loc[0]][tuple(g2loc[1:])]
        grid = refgenes.index(g)
        best_score = len(refgenes)
        best_loc = None
        for l in loc:
            lqid = non_overlap_id[l[0]][tuple(l[1:])]
            lqchr = l[0]
            score = abs(lqid - g1qid - grid + g1rid) + \
                abs(lqid - g2qid - grid + g2rid)
            if lqchr != g1loc[0] and lqchr != g2loc[0]:
                continue
            elif lqchr == g1loc[0] and lqchr != g2loc[0]:
                score += inter_chrom_penalty
            elif lqchr == g2loc[0] and lqchr != g1loc[0]:
                score += inter_chrom_penalty
            if score < best_score:
                best_score = score
                best_loc = l
        # Update
        if not best_loc:
            continue
        gene_syn[g] = best_loc
        gene_syn_assign[best_loc] = g
        related_genes = non_rb_assign[best_loc[0]][tuple(best_loc[1:])]
        for gene in related_genes:
            if gene[0] != g:
                if best_loc in non_rb_loc[gene[0]]:
                    non_rb_loc[gene[0]].remove(best_loc)
    syn_genes = list(gene_syn.keys())
    syn_genes.sort()
    with open('%s.gene.syn.tab' % prefix, 'w') as filep:
        filep.write('\n'.join(syn_genes))

    unassigned_genes = [
        g for g in refgenes if g not in syn_genes and g not in gene_rb]
    with open('%s.gene.unassign.tab' % prefix, 'w') as filep:
        filep.write('\n'.join(unassigned_genes))

    unassigned_seqs = [refseqs[refgenes.index(g)] for g in unassigned_genes]
    simple_fasta_write('%s.unassign.single.gene.fa' %
                       prefix, unassigned_genes, unassigned_seqs)

    # Save gff
    with open('%s.single.orf.gff' % prefix, 'w') as filep:
        content = ['##gff-version 3\n']
        for chrom in chroms:
            loc = list(non_overlap_id[chrom].keys())
            loc.sort()
            for i, l in enumerate(loc):
                start = l[0] + 1
                end = l[1]
                direct = '+' if l[2] == 1 else '-'
                note = 'ID=%s0%s%05dg' % (strain, chrom[3], i * 22 + 11)
                if l in rb_assign[chrom]:
                    note += ';Assign=%s;Method=ReciprocalBest' % rb_assign[chrom][l]
                elif (chrom, ) + l in gene_syn_assign:
                    note += ';Assign=%s;Method=Synteny' % gene_syn_assign[(
                        chrom, ) + l]
                else:
                    # Pick top 3 evidence
                    evi = sup_orf_loc[chrom][l][:3]
                    note += ';Assign=None;Hit='
                    hit = [e[0] for e in evi]
                    hit = ','.join(hit)
                    note += hit
                line = '%s\tCormackLab\tgene\t%d\t%d\t.\t%s\t.\t%s\n' % (
                    chrom, start, end, direct, note
                )
                content.append(line)
        filep.write(''.join(content))


def search_synteny(gene: str, genels: list)-> tuple:
    # genels.sort()
    p = 0
    gname_pattern = re.compile(r'[A-Z1-9]+0([A-Z])\d+g')
    gname = gname_pattern.match(gene).group(0)
    if gname < gname_pattern.match(min(genels)).group(0):
        return genels[0], genels[1], -1
    elif gname > gname_pattern.match(max(genels)).group(0):
        return genels[-2], genels[-1], 1
    for p in range(len(genels) - 1):
        g1 = genels[p]
        gn1 = gname_pattern.match(g1).group(0)
        g2 = genels[p + 1]
        gn2 = gname_pattern.match(g2).group(0)
        if gn1 < gene < gn2 or gn1 > gene > gn2:
            break
    gname_pattern = re.compile(r'[A-Z1-9]+0([A-Z])\d+g')
    chr1 = gname_pattern.match(g1).group(1)
    chr2 = gname_pattern.match(g2).group(1)
    chrq = gname_pattern.match(gene).group(1)
    if chr1 == chr2 == chrq:
        return g1, g2, 0
    elif chr1 == chrq and chrq != chr2:
        return genels[p-1], g1, 1
    elif chr2 == chrq and chrq != chr1:
        return g2, genels[p+1], -1
    else:
        print('Error Getting synteny info for %s, target %s, %s' % (gene, g1, g2))
        raise ValueError


def extract_single_cds(fa: str, gff: str, prefix: str):
    chroms, seqs = simple_fasta_load(fa)
    genome = {chroms[i]: seqs[i] for i in range(len(chroms))}
    genes = {}
    with open(gff) as filep:
        next(filep)
        for line in filep:
            ent = line.split()
            note = ent[8]
            note = note.split(';')
            id_ = note[0].split('=')[1]
            assign = note[1].split('=')[1]
            if assign == 'None':
                assign = note[2].split('=')[1]
                assign = assign.split(',')[0] + "_HIT"

            chrom = ent[0]
            start = int(ent[3]) - 1
            end = int(ent[4])
            direct = ent[6]
            cdna = genome[chrom][start:end]
            if direct == '-':
                cdna = rc_seq(cdna)
                gene = '%s/%s %d-%dC' % (id_, assign, end, start + 1)
            else:
                gene = '%s/%s %d-%dW' % (id_, assign, start, end)
            # Tmp: fix N->C
            if 'N' in cdna:
                cdna = re.sub('N', 'C', cdna)
                gene += ' N->C fix'
            cds = translate_exon(cdna)
            genes[gene] = [cdna, cds]

    gene_ls = list(genes.keys())
    gene_ls.sort()
    simple_fasta_write('%s.single.cdna.fa' %
                       prefix, gene_ls, [genes[g][0] for g in gene_ls])
    simple_fasta_write('%s.single.cds.fa' % prefix, gene_ls,
                       [genes[g][1] for g in gene_ls])


def extract_cdna(gff: str, fa: str, prefix: str):
    chroms, seqs = simple_fasta_load(fa)
    genome = {chroms[i]: seqs[i] for i in range(len(chroms))}
    gnames = []
    gseqs = []
    with open(gff) as filep:
        next(filep)
        for line in filep:
            ent = line.split()
            chrom = ent[0]
            start = int(ent[3])
            end = int(ent[4])
            direction = ent[6]
            seq = genome[chrom][start - 1: end]
            if direction == '-':
                seq = rc_seq(seq)
            note = ent[8]
            note = note.split(';')
            name = '%s/%s' % (note[0].split('=')[1], note[1].split('=')[1])
            gnames.append(name)
            gseqs.append(seq)
    simple_fasta_write(prefix + '.single.cdna.fa', gnames, gseqs)
