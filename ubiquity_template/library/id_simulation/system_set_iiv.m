function [cfg] = system_set_iiv(cfg, IIV1, IIV2, value)
% [cfg] = system_set_iiv(cfg, IIV1, IIV2, value)
%
% IIV1 and IIV2 name of the iiv block elements to set. If they are the same
% then the value represents the variance if they are different then the value
% represents the covariance.

isgood = true;
if(isfield(cfg, 'iiv'))

   IIV1_idx = strcmp(fieldnames(cfg.iiv.iivs), IIV1);
   IIV2_idx = strcmp(fieldnames(cfg.iiv.iivs), IIV2);
   if(sum(IIV1_idx) ~=1)
     vp(cfg,sprintf('IIV %s not found',IIV1));
     isgood = false;
   end
   if(sum(IIV2_idx) ~=1)
     vp(cfg,sprintf('IIV %s not found',IIV2));
     isgood = false;
   end

   if(isgood == true)
     cfg.iiv.values(IIV1_idx, IIV2_idx) = value;
     cfg.iiv.values(IIV2_idx, IIV1_idx) = value;
   else
     vp(cfg, 'system_set_iiv()');
     vp(cfg,'IIV information was not set');
     vp(cfg,'See above for more information');
   end

else
  vp(cfg, 'system_set_iiv()');
  vp(cfg, 'No IIV information was found') ;
  vp(cfg, 'These can be specified using:') ;
  vp(cfg, '<IIV:?>, <IIV:?:?>, and <IIVCOR:?:?>');
end
