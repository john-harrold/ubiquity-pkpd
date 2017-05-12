function [simtime] = system_ts_to_simtime(cfg, tstime, ts)


if(isfield(cfg.options.time_scales, ts))
  sim2tsscale = getfield(cfg.options.time_scales, ts);
   
  simtime = tstime./sim2tsscale;
 else
  vp(cfg, sprintf('Unable to find timescale %s', ts));
end

