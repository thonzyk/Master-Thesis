drop table if exists public.fus3_phospho_sites;

create table public.fus3_phospho_sites(
  id serial,
  yorf varchar(16),
  name varchar(16),
  position int,
  type varchar(1),
  primary key (id)
);