drop table if exists public.yorfs;

create table public.yorfs (
  id serial,
  yorf varchar(10),
  primary key (id)
);