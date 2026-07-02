import json

def main():
    with open('ssh-output/results/influenza-a-mp-indo-only/m05/meme_M1.json', 'r') as f:
        meme = json.load(f)
    meme_headers = meme['MLE']['headers']
    for idx, h in enumerate(meme_headers):
        print(f"Index {idx}: {h}")

if __name__ == '__main__':
    main()
