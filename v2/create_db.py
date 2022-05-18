import psycopg2
from psycopg2 import sql
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from getpass import getpass

with open('password.txt', 'r', encoding='utf-8') as fr:
    PSWD = fr.read()

# print('SQL server password: ')
# PSWD = getpass()

DB_NAME = 'masterthesis'


def a():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            pass


def create_the_db():
    con = psycopg2.connect(dbname='postgres', user='postgres', host='', password=PSWD)
    con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cur = con.cursor()
    cur.execute(f'drop database if exists {DB_NAME}')
    print(f'The {DB_NAME} database was dropped.')
    cur.execute(f"CREATE DATABASE {DB_NAME}")
    print(f'The {DB_NAME} database was created.')
    con.close()


if __name__ == '__main__':
    create_the_db()
