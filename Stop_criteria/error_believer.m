function obj=error_believer(search_x,candidate_x,kriging_model_lf,kriging_model_discrepancy,sample_xh,sample_xl,sample_yh,sample_yl,stop_value,cost_ratio)

% this function is used to determine the fidelity used to update the
% multi-fidelity surrogate model related to the estimated error
num_vari=size(search_x,2);
%%  if the candate sample is added in the LF sample set
fselected_lf=predictor(candidate_x,kriging_model_lf);
sample_xl_templf= [sample_xl;candidate_x];
sample_ly_templf = [sample_yl;fselected_lf];
% rebuild the Kriging model using the new design set
kriging_model_lf_temp=dacefit(sample_xl_templf,sample_ly_templf,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
sample_xh_yl_temp=predictor(sample_xh,kriging_model_lf_temp);
% rebuild the Kriging model using the new design set
kriging_model_discrepancy_lftemp=dacefit(sample_xh,sample_yh - sample_xh_yl_temp,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
% compute the entropy value
stop_value_lf=bootstrap_stop_multi_fidelity(search_x,kriging_model_lf_temp,kriging_model_discrepancy_lftemp,0.05);
% calculate the contribution of the lf model to the error reduction
lf_judge=cost_ratio*(stop_value-stop_value_lf)/stop_value;

%%  if the candate sample is added in the HF sample set
fselected_hf=estimated_value(candidate_x,kriging_model_lf,kriging_model_discrepancy);
sample_xh_temp = [sample_xh;candidate_x];
sample_hy_temp= [sample_yh;fselected_hf];
sample_ly_temphf=predictor(sample_xh_temp,kriging_model_lf);
% rebuild the Kriging model using the new design set
kriging_model_discrepancy_hftemp=dacefit(sample_xh_temp,sample_hy_temp - sample_ly_temphf,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
% compute the entropy value
stop_value_hf=bootstrap_stop_multi_fidelity(search_x,kriging_model_lf,kriging_model_discrepancy_hftemp,0.05);
% calculate the contribution of the lf model to the error reduction
hf_judge=(stop_value-stop_value_hf)/stop_value;

%% determine the fidelity used to add
if hf_judge<=0
  obj=1;
elseif   lf_judge>=hf_judge
     obj=1;
else
    obj=0;
end


end