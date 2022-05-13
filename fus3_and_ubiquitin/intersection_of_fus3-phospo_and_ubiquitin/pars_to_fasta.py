with open('data/ubi_fus3_pairs.txt', 'r', encoding='utf-8') as fr, open('data/ubi_fus3_pairs.fasta', 'w',
                                                                        encoding='utf-8') as fw:
    for line in fr:
        name, distance, fus3_pos, ubiq_pos, seq = line[:-1].split('\t')

        print(f'>{name}_distance={distance}_fus3_position={fus3_pos}_ubiq_position={ubiq_pos}', file=fw)
        print(seq, file=fw)
