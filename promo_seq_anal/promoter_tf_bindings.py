import numpy as np

import psycopg2
from promo_seq_anal.pfm_reader import read_pfm
from promo_seq_anal.dna_convolution import get_dna_convolution
from promo_seq_anal.constants import *
from matplotlib import pyplot as plt

CREATE_TABLE_TEC_SQL = '0010_create_tecbind_table.sql'
CREATE_TABLE_STE_SQL = '0011_create_stebind_table.sql'


def compute_responses(pfm_mat_list):
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            for pfm_file in pfm_mat_list:
                tf_name = pfm_file.split('.')[0]
                with open(RESULT_DIR / 'promoter_tf_response' / f'{tf_name}.tsv', 'w', encoding='utf-8') as fw:
                    print(f'gene_name\tmax_bind\tsum_bind\tresponse\trev_response', file=fw)
                    pfm = read_pfm(f'data/{pfm_file}', ignore_sum=True)
                    cur.execute('SELECT yorf, seq FROM public.promoters')

                    row = cur.fetchone()

                    while row is not None:
                        gene_name, dna = row
                        dna_conv = get_dna_convolution(dna, pfm)

                        with open(RESULT_DIR / 'promoter_tf_response' / tf_name / f'{gene_name}.tsv', 'w',
                                  encoding='utf-8') as fw2:
                            print(f'nucleotide\tresponse\tresponse_rev', file=fw2)
                            for ii in range(dna_conv.shape[1]):
                                nuc = dna[ii]
                                resp = dna_conv[0, ii]
                                resp_rev = dna_conv[1, ii]
                                print(f'{nuc}\t{resp}\t{resp_rev}', file=fw2)

                        max_binding = round(np.max(np.abs(dna_conv)), 4)
                        sum_binding = round(np.sum(np.abs(dna_conv)), 4)

                        response_str = ';'.join([str(el) for el in list(dna_conv[0, :])])
                        response_rev_str = ';'.join([str(el) for el in list(dna_conv[1, :])])

                        print(f'{gene_name}\t{max_binding}\t{sum_binding}\t{response_str}\t{response_rev_str}', file=fw)
                        row = cur.fetchone()


def run_responses_to_TCS_and():
    responses = compute_responses(['TEC1.pafm'])


def run_responses_to_STE12():
    responses = compute_responses(['STE12_consensus.pafm'])


def create_table_tecbind():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table...')
            cur.execute(open(f'sql/{CREATE_TABLE_TEC_SQL}', "r").read())
            with open(RESULT_DIR / 'promoter_tf_response' / 'TEC1.tsv', 'r', encoding='utf-8') as fr:
                fr.readline()
                print('Inserting data...')
                for i, line in enumerate(fr):
                    gene_name, max_bind, sum_bind, _, _ = line[:-1].split('\t')
                    max_bind = float(max_bind)
                    sum_bind = float(sum_bind)
                    assert sum_bind > max_bind

                    cur.execute("INSERT INTO public.tecbind (gene_name, max_bind, sum_bind) VALUES (%s, %s, %s)",
                                (gene_name, max_bind, sum_bind))

                    if i % 10000 == 0:
                        conn.commit()

    print('Creating initial indexes...')
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute("CREATE INDEX i_gene_name_tecbind ON public.tecbind (gene_name);")
            cur.execute("CREATE INDEX i_max_bind_tecbind ON public.tecbind (max_bind);")
            cur.execute("CREATE INDEX i_sum_bind_tecbind ON public.tecbind (sum_bind);")


def create_table_stebind():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table...')
            cur.execute(open(f'sql/{CREATE_TABLE_STE_SQL}', "r").read())
            with open(RESULT_DIR / 'promoter_tf_response' / 'STE12_consensus.tsv', 'r', encoding='utf-8') as fr:
                fr.readline()
                print('Inserting data...')
                for i, line in enumerate(fr):
                    gene_name, max_bind, sum_bind, _, _ = line[:-1].split('\t')
                    max_bind = float(max_bind)
                    sum_bind = float(sum_bind)
                    assert sum_bind >= max_bind

                    cur.execute("INSERT INTO public.stebind (gene_name, max_bind, sum_bind) VALUES (%s, %s, %s)",
                                (gene_name, max_bind, sum_bind))

                    if i % 10000 == 0:
                        conn.commit()

    print('Creating initial indexes...')
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute("CREATE INDEX i_gene_name_stebind ON public.stebind (gene_name);")
            cur.execute("CREATE INDEX i_max_bind_stebind ON public.stebind (max_bind);")
            cur.execute("CREATE INDEX i_sum_bind_stebind ON public.stebind (sum_bind);")


if __name__ == '__main__':
    # run_responses_to_STE12()
    create_table_stebind()
