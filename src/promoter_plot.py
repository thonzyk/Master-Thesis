import psycopg2
from matplotlib import pyplot as plt
from constants import *
import numpy as np


def plot_the_promoter(promoter_name):
    with psycopg2.connect(dbname=DB_NAME, user=USER, host="localhost", password=PSWD) as conn:
        with conn.cursor() as cur:
            cur.execute(
                'select pr.tec1_f, pr.tec1_r from public.promoter_response pr full join public.yorf_name yn on pr.yorf = yn.yorf  where yn.name = %s',
                (promoter_name,))
            result = cur.fetchall()
            if len(result) < 1:
                cur.execute(
                    'select pr.tec1_f, pr.tec1_r from public.promoter_response pr full join public.yorf_name yn on pr.yorf = yn.yorf  where yn.name = %s',
                    (promoter_name,))
                result = cur.fetchall()
            assert len(result) == 1
            result = result[0]
            resp_f = np.array(result[0])
            resp_r = -np.array(result[1])

            fig = plt.figure()
            plt.plot(resp_f, color='tab:blue')
            plt.plot(resp_r, color='tab:blue')
            plt.grid()
            # fig.axes[0].set_xticklabels([i for i in range(-1000, -8, 1)])
            plt.title(f'{promoter_name} promoter TEC1 binding')
            plt.xlabel('Position relative to start codon [bp]')
            plt.ylabel('Relative Binding Affinity')
            plt.show()


if __name__ == '__main__':
    plot_the_promoter('PRE10')
