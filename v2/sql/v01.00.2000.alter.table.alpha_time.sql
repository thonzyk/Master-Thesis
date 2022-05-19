alter table public.alpha_time drop column if exists avg_ratio;
alter table public.alpha_time add column if not exists avg_ratio float8;

alter table public.alpha_time drop column if exists avg_change;
alter table public.alpha_time add column if not exists avg_change float8;