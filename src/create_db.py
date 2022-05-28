import math

import pandas as pd
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from tqdm import tqdm

from constants import *
from utils import load_binding_mat, get_seq_response

# SQL scripts
CRATE_YORFS = 'sql/v01.00.0090.create.table.yorfs.sql'
CREATE_YORF_NAME = 'sql/v01.00.1000.create.table.yorf_name.sql'
CREATE_ALPHA_TIME = 'sql/v01.00.1010.create.table.alpha_time.sql'
CREATE_ALPHA_CONC = 'sql/v01.00.1020.create.table.alpha_conc.sql'
CREATE_REF_MRNA = 'sql/v01.00.1030.create.table.ref_mrna.sql'
CREATE_PROMOTERS = 'sql/v01.00.1040.create.table.promoters.sql'
CREATE_PROMOTER_RESPONSE = 'sql/v01.00.1050.create.table.promoter_response.sql'
ALTER_ALPHA_TIME = 'sql/v01.00.2000.alter.table.alpha_time.sql'


def add_columns():
    print('\n=== ADDING COLUMNS ===')

    # Alpha Time
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(ALTER_ALPHA_TIME, 'r', encoding='utf-8') as fr:
                # Add columns
                print('Adding columns to the table "alpha_time".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting data into the new columns of the "alpha_time" table.')
            cur.execute("select id, min_30, min_45, min_60, min_90, min_120 from public.alpha_time")

            result = cur.fetchall()

            for row in result:
                ID = row[0]
                values = []
                for item in row[1:]:
                    if item is not None:
                        values.append(item)

                assert len(values) > 1
                avg_resp = sum(values) / len(values)
                avg_change = sum([abs(values[i] - values[i + 1]) for i in range(len(values) - 1)]) / (len(values) - 1)
                cur.execute("update public.alpha_time set avg_ratio = %s, avg_change = %s where id=%s",
                            (avg_resp, avg_change, ID))

            print('Indexing the new columns of the "alpha_time" table.')
            for column in ['avg_ratio', 'avg_change']:
                cur.execute(f'create index i_alpha_time_{column} on public.alpha_time ({column})')

    # Promoters Bindings
    # Create "promoter_response" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_PROMOTER_RESPONSE, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "promoter_response".')
                cur.execute(fr.read())

    tec1_rbam = load_binding_mat()

    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            cur.execute("select yorf, seq from public.promoters")

            result = cur.fetchall()

            for row in tqdm(result):
                yorf = row[0]
                seq = row[1]

                promoter_response = get_seq_response(seq, tec1_rbam)

                cur.execute(
                    'insert into public.promoter_response (yorf, tec1_99, tec1_50, tec1_20, tec1_f, tec1_r) values (%s, %s, %s, %s, %s, %s)',
                    (yorf, promoter_response['score_99'], promoter_response['score_50'], promoter_response['score_20'],
                     promoter_response['response_f'], promoter_response['response_r']))

            # Index
            print('Indexing the table "promoter_response"')
            for column in ['tec1_99', 'tec1_50', 'tec1_20']:
                cur.execute(f'create index i_promoter_response_{column} on public.promoter_response ({column})')


def create_tables():
    print('\n=== TABLES CREATION ===')

    # Load Data
    all_yorfs = pd.read_csv('source_data/yeastract_com/yorfs.tsv', delimiter='\t')
    yorf_name_mapping = pd.read_csv('source_data/yeastract_com/yorf_name_mapping.tsv', delimiter='\t')
    alpha_time = pd.read_csv('source_data/Roberts_2000_PMID_10657304/2010.alpha_time.pcl', delimiter='\t')
    alpha_conc = pd.read_csv('source_data/Roberts_2000_PMID_10657304/2010.alpha_conc.pcl', delimiter='\t')
    ref_mrna = pd.read_csv('source_data/Lahtvee_2017_PMID_28365149/ref_mrna.tsv', delimiter='\t')
    promoters = pd.read_csv('source_data/ncbi/yeast_promoters.tsv', delimiter='\t')

    # Filter Data
    alpha_time = alpha_time[alpha_time['YORF'] != 'EWEIGHT']
    alpha_conc = alpha_conc[alpha_conc['YORF'] != 'EWEIGHT']

    # Create "yorfs" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CRATE_YORFS, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "yorfs".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "yorfs".')
            for i_row, row in all_yorfs.iterrows():
                yorf = row['yorf']
                cur.execute('insert into public.yorfs (yorf) values (%s)', (yorf,))

            # Index
            print('Indexing the table "yorfs"')
            for column in ['yorf']:
                cur.execute(f'create index i_yorfs_{column} on public.yorfs ({column})')

    # Create "yorf_name" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_YORF_NAME, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "yorf_name".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "yorf_name".')
            for i_row, row in yorf_name_mapping.iterrows():
                yorf = row['yorf']
                name = row['name']
                cur.execute('insert into public.yorf_name (yorf, name) values (%s, %s)', (yorf, name))

            # Index
            print('Indexing the table "yorf_name"')
            for column in ['yorf', 'name']:
                cur.execute(f'create index i_yorf_name_{column} on public.yorf_name ({column})')

    # Create "alpha_time" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_ALPHA_TIME, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "alpha_time".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "alpha_time".')
            for i_row, row in alpha_time.iterrows():
                yorf = row['YORF']
                min_0 = row['0 min'] if not math.isnan(row['0 min']) else None
                min_15 = row['15 min'] if not math.isnan(row['15 min']) else None
                min_30 = row['30 min'] if not math.isnan(row['30 min']) else None
                min_45 = row['45 min'] if not math.isnan(row['45 min']) else None
                min_60 = row['60 min'] if not math.isnan(row['60 min']) else None
                min_90 = row['90 min'] if not math.isnan(row['90 min']) else None
                min_120 = row['120 min'] if not math.isnan(row['120 min']) else None
                cur.execute(
                    'insert into public.alpha_time (yorf, min_0, min_15, min_30, min_45, min_60, min_90, min_120) values (%s, %s, %s, %s, %s, %s, %s, %s)',
                    (yorf, min_0, min_15, min_30, min_45, min_60, min_90, min_120))

            # Index
            print('Indexing the table "alpha_time"')
            for column in ['yorf', 'min_0', 'min_15', 'min_30', 'min_45', 'min_60', 'min_90', 'min_120']:
                cur.execute(f'create index i_alpha_time_{column} on public.alpha_time ({column})')

    # Create "alpha_conc" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_ALPHA_CONC, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "alpha_conc".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "alpha_conc".')
            for i_row, row in alpha_conc.iterrows():
                yorf = row['YORF']
                nm_0_15 = row['0.15 nM aF'] if not math.isnan(row['0.15 nM aF']) else None
                nm_0_5 = row['0.5 nM aF'] if not math.isnan(row['0.5 nM aF']) else None
                nm_1_5 = row['1.5 nM'] if not math.isnan(row['1.5 nM']) else None
                nm_5_0 = row['5 nM aF'] if not math.isnan(row['5 nM aF']) else None
                nm_15_8 = row['15.8 nM aF'] if not math.isnan(row['15.8 nM aF']) else None
                nm_50_0 = row['50 nM aF'] if not math.isnan(row['50 nM aF']) else None
                nm_158 = row['158 nM aF'] if not math.isnan(row['158 nM aF']) else None
                nm_500 = row['500 nM aF'] if not math.isnan(row['500 nM aF']) else None
                cur.execute(
                    'insert into public.alpha_conc (yorf, nm_0_15, nm_0_5, nm_1_5, nm_5_0, nm_15_8, nm_50_0, nm_158, nm_500) values (%s, %s, %s, %s, %s, %s, %s, %s, %s)',
                    (yorf, nm_0_15, nm_0_5, nm_1_5, nm_5_0, nm_15_8, nm_50_0, nm_158, nm_500))

            # Index
            print('Indexing the table "alpha_conc"')
            for column in ['yorf', 'nm_0_15', 'nm_0_5', 'nm_1_5', 'nm_5_0', 'nm_15_8', 'nm_50_0', 'nm_158',
                           'nm_500']:
                cur.execute(f'create index i_alpha_conc_{column} on public.alpha_conc ({column})')

    # Create "ref_mrna" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_REF_MRNA, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "ref_mrna".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "ref_mrna".')
            for i_row, row in ref_mrna.iterrows():
                yorf = row['YORF']
                mrna_in_pgdw = row['mRNA_in_pgDW']
                cur.execute('insert into public.ref_mrna (yorf, mrna_in_pgdw) values (%s, %s)',
                            (yorf, mrna_in_pgdw))

            # Index
            print('Indexing the table "ref_mrna"')
            for column in ['yorf', 'mrna_in_pgdw']:
                cur.execute(f'create index i_ref_mrna_{column} on public.ref_mrna ({column})')

    # Create "promoters" Table
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            with open(CREATE_PROMOTERS, 'r', encoding='utf-8') as fr:
                # Create
                print('Creating the table "promoters".')
                cur.execute(fr.read())

            # Insert Data
            print('Inserting source_data into the table "promoters".')
            for i_row, row in promoters.iterrows():
                yorf = row['yorf']
                seq = row['seq']
                cur.execute('insert into public.promoters (yorf, seq) values (%s, %s)', (yorf, seq))

            # Index
            print('Indexing the table "promoters"')
            for column in ['yorf', 'seq']:
                cur.execute(f'create index i_promoters_{column} on public.promoters ({column})')


def create_the_db():
    print('\n=== DATABASE CREATION ===')
    con = psycopg2.connect(dbname='postgres', user='postgres', host='', password=PSWD)
    con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cur = con.cursor()
    print(f'Dropping the "{DB_NAME}" database.')
    cur.execute(f'drop database if exists {DB_NAME}')
    print(f'Creating the "{DB_NAME}" database.')
    cur.execute(f"CREATE DATABASE {DB_NAME}")
    con.close()


if __name__ == '__main__':
    create_the_db()
    create_tables()
    add_columns()
