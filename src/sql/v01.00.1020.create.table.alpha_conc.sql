drop table if exists public.alpha_conc;

create table public.alpha_conc (
  id serial,
  yorf varchar(10),
  nm_0_15 float8,
  nm_0_5 float8,
  nm_1_5 float8,
  nm_5_0 float8,
  nm_15_8 float8,
  nm_50_0 float8,
  nm_158 float8,
  nm_500 float8,
  primary key (id)
);