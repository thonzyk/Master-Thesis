drop table if exists public.yeast_proteins;

create table public.yeast_proteins(
  id serial,
  yorf varchar(16),
  name varchar(16),
  seq varchar,
  primary key (id)
);