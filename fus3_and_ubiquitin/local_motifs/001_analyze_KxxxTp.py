import re

if __name__ == '__main__':
    with open('data/phospho_sites.txt', 'r', encoding='utf-8') as fr, open('results/K___T.txt', 'w',
                                                                           encoding='utf-8') as fw:
        for line in fr:
            line = line.replace('S@', 'S')
            if 'T@' not in line:
                continue

            res1 = re.search("K^[#]..T@", line)
            res2 = re.search("K#...T@", line)
            assert res1 is None or res2 is None

            if res1 is not None:
                result = res1.string[res1.regs[0][0]:res1.regs[0][1]]
            elif res2 is not None:
                result = res2.string[res2.regs[0][0]:res2.regs[0][1]]
            else:
                continue

            print(result.replace('#', ''), file=fw)

    with open('data/phospho_sites.txt', 'r', encoding='utf-8') as fr, open('results/K__T.txt', 'w',
                                                                           encoding='utf-8') as fw:
        for line in fr:
            line = line.replace('S@', 'S')
            if 'T@' not in line:
                continue

            res1 = re.search("K^[#].T@", line)
            res2 = re.search("K#..T@", line)
            assert res1 is None or res2 is None

            if res1 is not None:
                result = res1.string[res1.regs[0][0]:res1.regs[0][1]]
            elif res2 is not None:
                result = res2.string[res2.regs[0][0]:res2.regs[0][1]]
            else:
                continue

            print(result.replace('#', ''), file=fw)

    with open('data/phospho_sites.txt', 'r', encoding='utf-8') as fr, open('results/K____T.txt', 'w',
                                                                           encoding='utf-8') as fw:
        for line in fr:
            line = line.replace('S@', 'S')
            if 'T@' not in line:
                continue

            res1 = re.search("K^[#]...T@", line)
            res2 = re.search("K#....T@", line)
            assert res1 is None or res2 is None

            if res1 is not None:
                result = res1.string[res1.regs[0][0]:res1.regs[0][1]]
            elif res2 is not None:
                result = res2.string[res2.regs[0][0]:res2.regs[0][1]]
            else:
                continue

            print(result.replace('#', ''), file=fw)
