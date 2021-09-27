clc
clear all
close all
addpath('Test_function\test_function_multi_fidelity','Failure_Probability_evaluation','sample_plan_code','Infill_criteria','Stop_criteria')
%% defKine the tested function
%  multimodal_function  cubic_function  Circular_Pipe_function four_branches_function
for fun_name={'nonlinear_oscillator_function'}
    test_function=char(fun_name);
    % Cheap_multimodal_function  Cheap_cubic_function
    % Cheap_Cantilever_beam_function  Cheap_Circular_Pipe_function Cheap_four_branches_function
    test_lf_function=strcat('Cheap_',test_function);
    % test_lf_function='Cheap_multimodal_function';
    % import the basic information of the tested function
    [num_vari,z,A,cost_ratio,mu_vari,sigma_vari,design_space_vari,given_thresholds]=test_function_multi_fidelity(test_function);
    % define the basic parameters of the low fidelity function
    % the cost ratio
    mu=zeros(1,num_vari); sigma=ones(1,num_vari);
    design_space=[mu-5*sigma;mu+5*sigma];
    for error=given_thresholds
        num_trials=30;
        record.result=zeros(num_trials,8);
        for run=1:num_trials
            fprintf('--------------Test_function=%s;  run= %d  A=%d,  Cost_ratio=%d---------------------------\n',test_function,run,A,cost_ratio);
            % the paramters used for meta-model construction
            % LF samples
            num_low=4*cost_ratio;
            sample_xl= repmat(design_space(1,:),num_low,1)+repmat(design_space(2,:)-design_space(1,:),num_low,1).*...
                lhsdesign(num_low,num_vari,'criterion','maximin','iterations',1000);
            sample_yl= feval(test_lf_function,transfer_variables(sample_xl,mu_vari,sigma_vari),A);
            % HF samples
            num_high=6;
            [sample_xh,~]=subset(sample_xl,num_high);
            sample_yh = feval(test_function,transfer_variables(sample_xh,mu_vari,sigma_vari));
            %% constuct the surrogate model
            % low fidelity model
            kriging_model_lf=dacefit(sample_xl,sample_yl,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
            [sample_xh_yl,YYMSE1]=predictor(sample_xh,kriging_model_lf);
            % discrepance model
            kriging_model_discrepancy=dacefit(sample_xh,sample_yh - sample_xh_yl,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
            %% main iterative process
            gen=1;
            iteration=0;
            evaluation=size(sample_xh,1)+size(sample_xl,1)/cost_ratio;
            num_search=10^4*gen;
            search_x=MCS_Population_Generation(mu,sigma,num_search);
            [pf_estimate,pf_real,cov_estimate,real_relative_error]=reliabiliy_evaluation_multi_fidelity(search_x,mu,kriging_model_lf,kriging_model_discrepancy,test_function);
            while cov_estimate>0.05 ||iteration<=2
                num_search=10^4*gen;
                search_x=MCS_Population_Generation(mu,sigma,num_search);
                stop_value=bootstrap_stop_multi_fidelity(search_x,kriging_model_lf,kriging_model_discrepancy,0.05);
                while stop_value>error && iteration<500
                    %% step 1 determine the candidate sample point
                    best_mfru=Infill_MFRU(search_x,kriging_model_lf,kriging_model_discrepancy,sample_xh,sample_xl);
                    [bestobj,Index]=min(best_mfru);
                    candidate_x=search_x(Index,:);
                    if ismember(candidate_x,[sample_xl;sample_xh],'rows')==1
                        search_x=MCS_Population_Generation(mu,sigma,num_search);
                        best_mfru=Infill_MFRU(search_x,kriging_model_lf,kriging_model_discrepancy,sample_xh,sample_xl);
                        [bestobj,Index]=min(best_mfru);
                        candidate_x=search_x(Index,:);
                        fprintf('There is a sample coincide with the original ones!!! \n');
                    end
                    
                    %% step 2 determine the fidelity according to the error beleiver
                    fidelity=error_believer(search_x,candidate_x,kriging_model_lf,kriging_model_discrepancy,sample_xh,sample_xl,sample_yh,sample_yl,stop_value,cost_ratio);
                    
                    %% update the MF Kriging model
                    if fidelity==1
                        % add the new point to design set
                        if ismember(candidate_x,sample_xl,'rows')==0
                            sample_xl = [sample_xl;candidate_x];
                            Fselected_LF=feval(test_lf_function,transfer_variables(candidate_x,mu_vari,sigma_vari),A);
                            sample_yl = [sample_yl;Fselected_LF];
                            % rebuild the Kriging model using the new design set
                            kriging_model_lf=dacefit(sample_xl,sample_yl,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
                            sample_xh_yl=predictor(sample_xh,kriging_model_lf);
                            % rebuild the Kriging model using the new design set
                            kriging_model_discrepancy=dacefit(sample_xh,sample_yh - sample_xh_yl,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
                        else
                            search_x=MCS_Population_Generation(mu,sigma,num_search);
                        end
                    elseif fidelity==0
                        % add the new point to design set
                        if ismember(candidate_x,sample_xh,'rows')==0
                            sample_xh = [sample_xh;candidate_x];
                            Fselected_HF=feval(test_function,transfer_variables(candidate_x,mu_vari,sigma_vari));
                            [Fselected_LF,~]=predictor(candidate_x,kriging_model_lf);
                            sample_yh = [sample_yh;Fselected_HF];
                            sample_xh_yl = [sample_xh_yl;Fselected_LF];
                            % rebuild the Kriging model using the new design set
                            kriging_model_discrepancy=dacefit(sample_xh,sample_yh - sample_xh_yl,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
                            
                        else
                            search_x=MCS_Population_Generation(mu,sigma,num_search);
                        end
                    else
                        fprintf('Something wrong!!! \n');
                    end
                    % updating some parameters
                    iteration = iteration+1;
                    %                     if size(sample_xh,2)==2
                    %                         plot_figures_multi_fidelity(kriging_model_lf, kriging_model_discrepancy,design_space,num_low,num_high,sample_xl,sample_xh,test_hf_function,iteration);
                    %                     end
                    evaluation = size(sample_xl,1)/cost_ratio+size(sample_xh,1);
                    %%
                    %the stopping criterion
                    stop_value=bootstrap_stop_multi_fidelity(search_x,kriging_model_lf,kriging_model_discrepancy,0.05);
                    [pf_estimate,pf_real,~,real_relative_error]=reliabiliy_evaluation_multi_fidelity(search_x,mu,kriging_model_lf,kriging_model_discrepancy,test_function);
                    fprintf(' A=%f,cr=%d,Run=%d; iter=%d; tc_hf=%d,tc_lf=%d, pf_estimate=%f; pf_real=%f; stop_value=%f  real_error=%f \n',...
                        A,cost_ratio,run,iteration ,size(sample_xh,1), size(sample_xl,1),pf_estimate,pf_real,stop_value,real_relative_error);
                end
                %% calclulate the failure probility
                [pf_estimate,pf_real,cov_estimate,real_relative_error]=reliabiliy_evaluation_multi_fidelity(search_x,mu,kriging_model_lf,kriging_model_discrepancy,test_function);
                % final result of one run
                gen=gen+1;
                high_size=size(sample_xh,1);
                low_size=size(sample_xl,1);
                cost=high_size+low_size/cost_ratio;
                
            end
            record.result(run,:)=[high_size,low_size,cost,pf_estimate,pf_real,cov_estimate,real_relative_error,size(search_x,1)];
            fprintf('A=%f,cr=%d,Run=%d; iter=%d; tc_hf=%d,tc_lf=%d, pf_estimate=%f; pf_real=%f; stop_value=%f  real_error=%f \n',...
                A,cost_ratio,run,iteration , size(sample_xh,1),  size(sample_xl,1),pf_estimate,pf_real,stop_value,real_relative_error);
            save(strcat('Results/',test_function,'/',mfilename,'_',num2str(100*error),'_',num2str(cost_ratio),'.mat'),'record');
        end
        record.mean=mean(record.result);
        save(strcat('Results/',test_function,'/',mfilename,'_',num2str(100*error),'_',num2str(cost_ratio),'.mat'),'record');
        
    end
    
    
end

