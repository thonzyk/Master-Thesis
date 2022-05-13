import pandas as pd
import psycopg2
import numpy as np


def find_closest_pairs():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute('SELECT yorf, name, position, type FROM public.fus3_phospho_sites')
            result = cur.fetchall()

            for row in result:
                yorf, name, position, phospho_type = row
                position_0 = position - 1

                cur.execute('SELECT position FROM public.ubiquitin_sites where yorf=%s', (yorf,))
                result = cur.fetchall()

                if len(result) < 1:
                    continue

                result_positions = [item[0] for item in result]
                distances = [abs(result_positions[i] - position) for i in range(len(result_positions))]

                closest_match = np.argmin(distances)
                closest_position = result_positions[closest_match]
                closest_distance = min(distances)

                cur.execute(
                    'INSERT INTO public.phospho_ubiq_pairs (yorf, name, position_phospho, position_ubiq, distance, source) values (%s, %s, %s, %s, %s, %s)',
                    (yorf, name, position, closest_position, closest_distance, phospho_type))


if __name__ == '__main__':
    find_closest_pairs()
