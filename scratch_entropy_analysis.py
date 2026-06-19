import yaml

with open('ssh-output/influenza-a-mp-universal-retake/m07/primer_targets_manifest.yaml') as f:
    manifest = yaml.safe_load(f)

region = manifest['candidate_regions'][0]
entropy = region['positional_scores']['entropy']
meme = region['positional_scores']['meme_pval']
total = len(entropy)

below005 = sum(1 for x in entropy if x < 0.05)
below010 = sum(1 for x in entropy if x < 0.10)
below020 = sum(1 for x in entropy if x < 0.20)
above020 = sum(1 for x in entropy if x >= 0.20)
above010 = sum(1 for x in entropy if x >= 0.10)

print(f"Total positions: {total}")
print(f"< 0.05 (sangat konservatif): {below005} ({below005/total*100:.1f}%)")
print(f"< 0.10 (konservatif)       : {below010} ({below010/total*100:.1f}%)")
print(f"< 0.20 (aman untuk primer) : {below020} ({below020/total*100:.1f}%)")
print(f">= 0.20 (akan dikecualikan): {above020} ({above020/total*100:.1f}%)")
print(f">= 0.10 (dikecualikan MEME): {above010} ({above010/total*100:.1f}%)")
print()

# Cari blok aman (entropy < 0.20, panjang >= 18bp untuk primer min length)
MIN_PRIMER_LEN = 18
safe_blocks = []
in_block = False
block_start = 0
for i, val in enumerate(entropy):
    if val < 0.20:
        if not in_block:
            in_block = True
            block_start = i
    else:
        if in_block:
            length = i - block_start
            if length >= MIN_PRIMER_LEN:
                safe_blocks.append((block_start, i-1, length))
            in_block = False
if in_block:
    length = total - block_start
    if length >= MIN_PRIMER_LEN:
        safe_blocks.append((block_start, total-1, length))

print(f"Safe blocks (entropy < 0.20, len >= {MIN_PRIMER_LEN}bp):")
for s, e, l in safe_blocks:
    avg_e = sum(entropy[s:e+1]) / l
    print(f"  pos {s:4d}-{e:4d} | len={l:4d}bp | avg_entropy={avg_e:.4f}")

print()
print(f"Total safe space: {sum(l for _,_,l in safe_blocks)}bp out of {total}bp")
print(f"Number of safe blocks: {len(safe_blocks)}")
