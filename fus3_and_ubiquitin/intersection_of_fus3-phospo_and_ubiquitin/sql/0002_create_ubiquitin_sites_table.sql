drop table if exists public.ubiquitin_sites;

create table public.ubiquitin_sites(
  id serial,
  yorf varchar(16),
  name varchar(16),
  position int,
  source varchar(32),
  primary key (id)
);