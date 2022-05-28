drop table if exists public.promoters;

create table public.promoters(
  id serial,
  yorf varchar(10),
  seq varchar(1000),
  primary key (id)
);