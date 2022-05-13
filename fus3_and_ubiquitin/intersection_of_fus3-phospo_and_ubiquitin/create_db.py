import pandas as pd
import psycopg2

CREATE_FUS3_PHOSPHO_SITES_TABLE_SQL = '0001_create_fus3_phospho_sites_table.sql'
CREATE_UBIQUITIN_SITES_TABLE_SQL = '0002_create_ubiquitin_sites_table.sql'

CREATE_PHOSPHO_UBIQ_PAIRS_TABLE_SQL = '0003_create_phospho_ubiq_pairs_table.sql'
CREATE_YEAST_PROTEINS_TABLE_SQL = '0004_create_yeast_proteins_table.sql'


def create_tables():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table fus3_phospho_sites...')
            cur.execute(open(f'sql/{CREATE_FUS3_PHOSPHO_SITES_TABLE_SQL}', "r").read())

    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table ubiquitin_sites...')
            cur.execute(open(f'sql/{CREATE_UBIQUITIN_SITES_TABLE_SQL}', "r").read())

    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table phospho_ubiq_pairs...')
            cur.execute(open(f'sql/{CREATE_PHOSPHO_UBIQ_PAIRS_TABLE_SQL}', "r").read())

    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating table yeast_proteins...')
            cur.execute(open(f'sql/{CREATE_YEAST_PROTEINS_TABLE_SQL}', "r").read())

    print()


def index_tables():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            print('Creating initial indexes for fus3_phospho_sites...')
            cur.execute("CREATE INDEX i_fus3_phospho_sites_yorf ON public.fus3_phospho_sites (yorf);")
            cur.execute("CREATE INDEX i_fus3_phospho_sites_name ON public.fus3_phospho_sites (name);")
            cur.execute("CREATE INDEX i_fus3_phospho_sites_position ON public.fus3_phospho_sites (position);")

            print('Creating initial indexes for ubiquitin_sites...')
            cur.execute("CREATE INDEX i_ubiquitin_sites_yorf ON public.ubiquitin_sites (yorf);")
            cur.execute("CREATE INDEX i_ubiquitin_sites_name ON public.ubiquitin_sites (name);")
            cur.execute("CREATE INDEX i_ubiquitin_sites_position ON public.ubiquitin_sites (position);")
            cur.execute("CREATE INDEX i_ubiquitin_sites_source ON public.ubiquitin_sites (source);")

            print('Creating initial indexes for phospho_ubiq_pairs...')
            cur.execute("CREATE INDEX i_phospho_ubiq_pairs_yorf ON public.phospho_ubiq_pairs (yorf);")
            cur.execute("CREATE INDEX i_phospho_ubiq_pairs_name ON public.phospho_ubiq_pairs (name);")
            cur.execute(
                "CREATE INDEX i_phospho_ubiq_pairs_position_phospho ON public.phospho_ubiq_pairs (position_phospho);")
            cur.execute("CREATE INDEX i_phospho_ubiq_pairs_position_ubiq ON public.phospho_ubiq_pairs (position_ubiq);")
            cur.execute("CREATE INDEX i_phospho_ubiq_pairs_source ON public.phospho_ubiq_pairs (source);")

            print('Creating initial indexes for yeast_proteins...')
            cur.execute("CREATE INDEX i_yeast_proteins_yorf ON public.yeast_proteins (yorf);")
            cur.execute("CREATE INDEX i_yeast_proteins_name ON public.yeast_proteins (name);")

    print()


def cz_str_to_float(item):
    return float(item.replace(',', '.'))


def insert_into_tables():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            data = pd.read_csv('data/S288C_proteins.tsv', delimiter='\t', index_col=False)
            print('Inserting protein sequences into yeast_proteins...')
            for row_i, row in data.iterrows():
                yorf, name, seq = row['locus_name'], row['gene_name'], row['seq']

                cur.execute("INSERT INTO public.yeast_proteins (yorf, name, seq) VALUES (%s, %s, %s)",
                            (yorf, name, seq))

    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            data = pd.read_csv('data/fus3_phospho/2011_Ronghu.tsv', delimiter='\t', index_col=False)
            print('Inserting protein sequences into fus3_phospho_sites...')
            for row_i, row in data.iterrows():

                yorf, name = row['Reference'], row['Gene symbol']
                positions = [row['Site 1'], row['Site 2'], row['Site 3']]
                a_scores = [cz_str_to_float(row['Ascore 1']), cz_str_to_float(row['Ascore 2']),
                            cz_str_to_float(row['Ascore 3'])]

                for i in range(len(positions)):
                    if positions[i] < 1:
                        continue

                    if a_scores[i] < 13.0:
                        continue

                    cur.execute("INSERT INTO public.fus3_phospho_sites (yorf, name, position ) VALUES (%s, %s, %s)",
                                (yorf, name, positions[i]))

    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            data = pd.read_csv('data/ubiquitin/2013_Swaney.tsv', delimiter='\t', index_col=False)
            print('Inserting protein sequences into ubiquitin_sites...')
            for row_i, row in data.iterrows():

                yorf, name, reside, position, a_score = row['Protein'], None, row['Reside'], row['Position'], row[
                    'Ascore']

                a_score = str(a_score).replace(',', '.')
                if not a_score.isnumeric():
                    continue
                a_score = float(a_score)

                if a_score < 13.0:
                    continue

                cur.execute(
                    "INSERT INTO public.ubiquitin_sites (yorf, name, position, source ) VALUES (%s, %s, %s, %s)",
                    (yorf, name, position, '2013_Swaney'))


def add_phospho_type():
    n_errors = 0
    print('Adding phosphorylated residues...')
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute('SELECT yorf, name, position FROM public.fus3_phospho_sites')
            result = cur.fetchall()

            for ros_i, row in enumerate(result):
                yorf, name, position = row
                position_0 = position - 1  # 1-indexing to 0-indexing

                cur.execute('SELECT seq FROM public.yeast_proteins where yorf=%s', (yorf,))

                seqs = cur.fetchall()
                assert len(seqs) == 1
                seq = seqs[0][0]

                if position_0 > len(seq) - 1:
                    n_errors += 1
                    continue

                phospho_residue = seq[position_0]

                if not (phospho_residue == 'T' or phospho_residue == 'S' or phospho_residue == 'Y'):
                    print(f'unexpected phosphorylated amino acid: {phospho_residue}')

                cur.execute('UPDATE public.fus3_phospho_sites SET type=%s where yorf=%s', (phospho_residue, yorf))

    print(f'Errors: {n_errors}')


if __name__ == '__main__':
    create_tables()
    index_tables()
    insert_into_tables()
    add_phospho_type()
