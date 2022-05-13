drop table if exists public.phospho_ubiq_pairs;

create table public.phospho_ubiq_pairs(
  id serial,
  yorf varchar(16),
  name varchar(16),
  position_phospho int,
  position_ubiq int,
  distance int,
  source varchar(32),
  primary key (id)
);