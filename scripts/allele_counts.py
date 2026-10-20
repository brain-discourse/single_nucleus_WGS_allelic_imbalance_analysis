#!/usr/bin/env python3
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-18h-54m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: python script for mpileup parsing and counting
# ==========================================================
import re
import sys

comp = re.compile(r'[\^].')  # start of read + mapq


def clean(b: str) -> str:
    b = comp.sub('', b)
    b = b.replace('$', '')
    i = 0
    out = []
    while i < len(b):
        if b[i] in '+-':
            j = i + 1
            while j < len(b) and b[j].isdigit():
                j += 1
            if j == i + 1:
                i += 1
                continue
            n = int(b[i + 1:j])
            i = j + n
        else:
            out.append(b[i])
            i += 1
    return ''.join(out)


def main():
    print("\t".join(["chrom", "pos", "ref", "dp", "A", "C", "G", "T"]))

    for line in sys.stdin:
        if not line.strip():
            continue
            
        chrom, pos, ref, dp_str, bases, _ = line.rstrip('\n').split('\t')[:6]
        dp = int(dp_str)
        b = clean(bases)
        
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        ref_upper = ref.upper()
        
        for ch in b:
            if ch in '.,':
                if ref_upper in counts:
                    counts[ref_upper] += 1
            else:
                ch_upper = ch.upper()
                if ch_upper in counts:
                    counts[ch_upper] += 1
                    
        sys.stdout.write("\t".join([
            chrom, pos, ref, str(dp),
            str(counts['A']), str(counts['C']),
            str(counts['G']), str(counts['T'])
        ]) + "\n")


if __name__ == "__main__":
    main()
