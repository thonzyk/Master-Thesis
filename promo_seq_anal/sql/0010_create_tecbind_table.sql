drop table if exists public.tecbind;

create table public.tecbind(
  id serial,
  gene_name varchar,
  max_bind float,
  sum_bind float,
  primary key (id)
);