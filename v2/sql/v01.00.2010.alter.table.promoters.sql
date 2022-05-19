alter table public.promoters drop column if exists tec1_99;
alter table public.promoters add column if not exists tec1_99 int;

alter table public.promoters drop column if exists tec1_50;
alter table public.promoters add column if not exists tec1_50 int;

alter table public.promoters drop column if exists tec1_25;
alter table public.promoters add column if not exists tec1_25 int;