import psycopg2


def add_the_avg_column_in_time_response():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute("select id, min_30, min_45, min_60, min_90, min_120 from public.alpha_factor_time_response")

            result = cur.fetchall()

            for row in result:
                ID = row[0]
                cnt = 0
                sm = 0.0
                for item in row[1:]:
                    if item is not None:
                        cnt += 1
                        sm += item

                avg_resp = sm / cnt
                cur.execute("update public.alpha_factor_time_response set avg_rate = %s where id=%s", (avg_resp, ID))


def add_the_avg_column_in_conc_response():
    with psycopg2.connect(dbname="masterthesis", user="postgres", host="localhost", password="mlBIO") as conn:
        with conn.cursor() as cur:
            cur.execute(
                "select id, af_0_15_nm, af_0_5_nm, af_1_5_nm, af_5_0_nm, af_15_8_nm, af_50_0_nm, af_158_nm, af_500_nm from public.alpha_factor_concentration_response")

            result = cur.fetchall()

            for row in result:
                ID = row[0]
                cnt = 0
                sm = 0.0
                for item in row[1:]:
                    if item is not None:
                        cnt += 1
                        sm += item

                avg_resp = sm / cnt
                cur.execute("update public.alpha_factor_concentration_response set avg_rate = %s where id=%s",
                            (avg_resp, ID))


if __name__ == '__main__':
    add_the_avg_column_in_conc_response()
