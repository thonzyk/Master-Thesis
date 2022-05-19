drop table if exists public.yorf_name;

create table public.yorf_name (
  id serial,
  yorf varchar(10),
  name varchar(50),
  primary key (id)
);