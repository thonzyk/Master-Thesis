USER = 'postgres'
with open('password.txt', 'r', encoding='utf-8') as fr:
    PSWD = fr.read()

# print('SQL server password: ')
# PSWD = getpass()

DB_NAME = 'masterthesis'