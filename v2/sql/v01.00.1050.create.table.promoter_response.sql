drop table if exists public.promoter_response;

create table public.promoter_response(
  id serial,
  yorf varchar(10),
  tec1_99 int4,
  tec1_50 int4,
  tec1_20 int4,
  tec1_f float4 array[992],
  tec1_r float4 array[992],
  primary key (id)
);