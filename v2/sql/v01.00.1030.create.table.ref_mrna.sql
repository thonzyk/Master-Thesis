drop table if exists public.ref_mrna;

create table public.ref_mrna (
  id serial,
  yorf varchar(10),
  mrna_in_pgdw float8,
  primary key (id)
);