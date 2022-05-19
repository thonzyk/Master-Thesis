drop table if exists public.alpha_time;

create table public.alpha_time (
  id serial,
  yorf varchar(10),
  min_0 float8,
  min_15 float8,
  min_30 float8,
  min_45 float8,
  min_60 float8,
  min_90 float8,
  min_120 float8,
  primary key (id)
);