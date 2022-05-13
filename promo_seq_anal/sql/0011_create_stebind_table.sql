drop table if exists public.stebind;

create table public.stebind(
  id serial,
  gene_name varchar,
  max_bind float,
  sum_bind float,
  primary key (id)
);