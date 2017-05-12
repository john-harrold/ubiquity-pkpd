function [s]=solution_statistics(parameters, cfg)
% function [s]=solution_statistics(parameters, cfg)
%
% Function to calculate the solution statistics or a set of 'parameters'. The
% configuration variable 'cfg' is a data structure. See the help for
% 'check_cfg' for a listing of the required fields
%
% The solution statistics are stored in the structure 's'. It has the
% following fields:
%
% General fields:
%   s.coefficient_of_variation                CV percent
%   s.confidence_interval.lower_bound         Lower bound on 95% confidence interval
%   s.confidence_interval.upper_bound         Upper bound on 95% confidence interval
%   s.covariance                              Covariance matrix
%   s.degrees_of_freedom                      Number of degrees of freedom
%   s.num_observations                        Total number of observations
%   s.aic                                     Akiake Information Criterion
%   s.bic                                     Bayesian Information Criterion
%
% Intermediate values used and returned 
% For the 'wls' analysis type
%   s.wls.jacobian                            Full Jacobian matrix
%   s.wls.weights                             Weighting matrix
%   s.wls.error_variance                      Error variance matrix
%
% For the 'ml' analysis type
%   s.ml.M                                    M matrix 
%
% Changelog:
%    2010.12.15     decreased the values used for relative and 
%                   absolute tolerance by two orders of magnitude each
%                   to improve solution statistics calculation



% perturbing the parameters for numerical calculation of derivatives
% RelTol and Abstol are the defaults returned by simset
%RelTol = 1e-3;
%AbsTol = 1e-6;
RelTol = 1e-5;
AbsTol = 1e-8;
perturbation = max(abs(parameters), AbsTol).*RelTol;

% 
% initializing variables
% 

if(strcmp('wls', cfg.estimation.objective_type))
    s.wls       = {};
elseif(strcmp('ml', cfg.estimation.objective_type))
    s.ml        = {};
end



 try 
  % obtaining observation informaiotn by
  % evaluating the user defined function
  eval(sprintf('[observations]=%s(parameters,cfg);',cfg.estimation.observation_function));

  %
  % calculating the objective function value
  %
  [objective]=calculate_objective(parameters,cfg);

  %
  % perturbing the parameters for calculating the
  % partial derivatives
  %

  for i=1:length(parameters);
    % clearing the contents of perturb_tmp
    clear perturb_tmp;
    perturb_tmp.parameters    = parameters;
    perturb_tmp.parameters(i) = perturb_tmp.parameters(i) +  perturbation(i) ;

    % finding the perturbed observations
    eval(sprintf('[perturb_tmp.observations]=%s(perturb_tmp.parameters,cfg);',cfg.estimation.observation_function));
  
    % saving the contents of perturb_tmp into the 
    % perturbations data structure
    eval(sprintf('perturbations_plus.parameter%d = perturb_tmp;',i));

    clear perturb_tmp;
    perturb_tmp.parameters    = parameters;
    perturb_tmp.parameters(i) = perturb_tmp.parameters(i) -  perturbation(i) ;

    % finding the perturbed observations
    eval(sprintf('[perturb_tmp.observations]=%s(perturb_tmp.parameters,cfg);',cfg.estimation.observation_function));

    % saving the contents of perturb_tmp into the 
    % perturbations data structure
    eval(sprintf('perturbations_minus.parameter%d = perturb_tmp;',i));

  end
  catch
    disp(sprintf(' -> Unable to retrieve observations'));
    disp(sprintf(' -> This may result if the system is stiff around the parameter estimate.'));
    return
  end


  % notes:
  % index lengths:
  %  p - number of parameters
  %  q - number of variance parameters
  %  l - number of outputs
  %  m - number of time points 

  % getting all subjects
  all_subjects          = fieldnames(observations);
  all_outputs           = {};
  num_measurements.all  = 0; % total number of measurements

  % 
  % analyzing 'observations' to determine the number of possible outputs, time
  % points etc.
  % 
  for  subidx =1:length(all_subjects)
    % pulling the individual assocaited with subid
    individual = getfield(observations, char(all_subjects(subidx)));

    % getting that individuals outputs
    outputs  = fieldnames(individual);
    % now looping through each output in that subject
    for  outidx =1:length(outputs)
      if(sum(strcmp(all_outputs, char(outputs(outidx)))) < 1)
          all_outputs = [all_outputs; char(outputs(outidx))];
          eval(sprintf('num_measurements.outputs.%s = 0;', char(outputs(outidx)) ));
      end
      eval(sprintf('[tmpsize] = size(observations.%s.%s);', ...
                   char(all_subjects(subidx)), ...
                   char(outputs(outidx))));
      % calculating the number of sample times
      num_measurements.all = num_measurements.all  +  tmpsize(1);
      eval(sprintf('num_measurements.outputs.%s = num_measurements.outputs.%s + tmpsize(1) ;', ... 
          char(outputs(outidx)), ...
          char(outputs(outidx)) ));
    end
  end

  %sorting the lists of subjects and outputs
  all_subjects = sort(all_subjects);
  all_outputs  = sort(all_outputs);

  % calcuating the solution statistics for the 
  % weight least squares option
  if(strcmp('wls', cfg.estimation.objective_type))

    % temporary variables to contain the Jacobian (P), matrix of weights (W),
    % and the error variance matrix (G)
    P            = zeros(num_measurements.all, length(parameters));
    W            = zeros(num_measurements.all, 1);
    G            = [];
    time_counter = 0;

    for outidx = 1:length(all_outputs)
      % all predictions, obsevations, and 
      % weights for this output
      pred_output = [];
      obs_output  = [];
      wt_output   = [];
      eval(sprintf('num_output   = num_measurements.outputs.%s;',  char(all_outputs(outidx))));
      df_output   = num_output - (length(parameters)./length(all_outputs));

      for subidx = 1:length(all_subjects)
         % fetching data for the subject
         subject = getfield(observations,  char(all_subjects(subidx)));

         % checking to see if the subject has the current output
         if(isfield(subject, char(all_outputs(outidx))))
           output         = getfield(subject,  char(all_outputs(outidx)));
           tmpsize        = size(output);
           time_offset    = tmpsize(1);
           pred_output    = [pred_output; output(:,3)];
           obs_output     = [obs_output;  output(:,2)];
           wt_output      = [wt_output;   output(:,4)];

           W((time_counter+1):(time_counter+time_offset),1) = 1./output(:,4); 
           
           for pidx = 1:length(parameters)
            eval(sprintf('perturb_tmp = perturbations_plus.parameter%d;',pidx));

            deltaoutput_plus  =    getfield(getfield(perturb_tmp.observations,...
                                          char(all_subjects(subidx))), ...
                                          char(all_outputs(outidx)));

            eval(sprintf('perturb_tmp = perturbations_minus.parameter%d;',pidx));

            deltaoutput_minus =    getfield(getfield(perturb_tmp.observations,...
                                                char(all_subjects(subidx))), ...
                                                char(all_outputs(outidx)));
            
            % calculating the partial derivative 
            partiald = (deltaoutput_plus(:,3) - deltaoutput_minus(:,3))./(2*perturbation(pidx));

            P((time_counter+1):(time_counter+time_offset),pidx) = partiald;
          
           end
           time_counter = time_counter + time_offset;
         end
      end
      % calculating the variance for the current output 
      % and making an entry in the G vector for these values
      variance_output = sum(1./wt_output.*((pred_output-obs_output).^2))./df_output;
      G = [G; ones(num_output,1).*variance_output];

    end

    % Creating a diagional matrix from W and G
    % For example:
    %
    %                        | 1     0     0  |
    % W = [1 2 3]  ----> W = | 0     2     0  |
    %                        | 0     0     3  |
     
    W = diag(W);
    G = diag(G);


    s.wls.jacobian           = P;
    s.wls.weights            = W;
    s.wls.error_variance     = G;
    s.covariance             = inv(P'*W*P)*(P'*W*G*W*P)*inv(P'*W*P);
    s.degrees_of_freedom     = num_measurements.all - length(parameters);
    s.num_observations       = num_measurements.all;
    s.objective              = objective;
    s.aic                    = s.num_observations*log(objective) + 2*length(parameters);
    s.bic                    = s.num_observations*log(objective) + log(s.num_observations)*length(parameters);

  % calcuating the solution statistics 
  % for maximum likelihood option
  elseif(strcmp('ml', cfg.estimation.objective_type))
    % M has the following structure:
    %     _                                 _    
    %    |                     ^             |        
    %    |                     :             |        
    %    |                     :             |        
    %    |         MI        p :     MIII    |        
    %    |                     :             |        
    %    |                     :             |        
    %    |       p             :             |        
    %    |<.................................>|        
    %    |                     :      q      |        
    %    |                     :             |        
    %    |         MIII        :q    MII     |        
    %    |                     :             |        
    %    |_                    v            _|        
    %                                                     
    %    Where p is the number of system parameters and 
    %    q is the number of variance parameters
    %                                                     
    %    With the three block components MI, MII and MIII     

    % initializing M
    M     = zeros(length(parameters), length(parameters));
    dim.p = cfg.estimation.parameters.system;
    dim.q = length(parameters) - dim.p;

    for j=1:length(parameters)
      for k=1:length(parameters)
        %looping through each subject and output
        for outidx = 1:length(all_outputs)
          for subidx = 1:length(all_subjects)
            % fetching data for the subject
            subject = getfield(observations,  char(all_subjects(subidx)));
            % checking to see if the subject has the current output
            if(isfield(subject, char(all_outputs(outidx))))
              % clearing temporary variables used at each pass through
              clear partials perturbj_plus perturbk_plus perturbj_minus perturbk_minus;
              clear outputj outputk output outg;

              eval(sprintf('perturbj_plus = perturbations_plus.parameter%d;',j));
              eval(sprintf('perturbk_plus = perturbations_plus.parameter%d;',k));
              eval(sprintf('perturbj_minus= perturbations_minus.parameter%d;',j));
              eval(sprintf('perturbk_minus= perturbations_minus.parameter%d;',k));
              % outputs at parameter estimates
              output      = getfield(subject,  char(all_outputs(outidx)));
              outg        = output(:,4);

              outputj_plus= getfield(getfield(perturbj_plus.observations,...
                                     char(all_subjects(subidx))), ...
                                     char(all_outputs(outidx)));
              outputk_plus= getfield(getfield(perturbk_plus.observations,...
                                         char(all_subjects(subidx))), ...
                                         char(all_outputs(outidx)));

              outputj_minus= getfield(getfield(perturbj_minus.observations,...
                                     char(all_subjects(subidx))), ...
                                     char(all_outputs(outidx)));
              outputk_minus= getfield(getfield(perturbk_minus.observations,...
                                         char(all_subjects(subidx))), ...
                                         char(all_outputs(outidx)));
               
              partials.gj = (outputj_plus(:,4)-outputj_minus(:,4))./(2*perturbation(j));
              partials.gk = (outputk_plus(:,4)-outputk_minus(:,4))./(2*perturbation(k));
              

              partials.yj = (outputj_plus(:,3)-outputj_minus(:,3))./(2*perturbation(j));
              partials.yk = (outputk_plus(:,3)-outputk_minus(:,3))./(2*perturbation(k));

              if (j <= dim.p)& (k <=dim.p)
                %                                                     
                % Section MI
                %                                                     

                M(j,k) = M(j,k) + 1/2*sum(partials.gj.*partials.gk./outg.^2);
                M(j,k) = M(j,k) +     sum(partials.yj.*partials.yk./outg);
              else
                %                                                     
                % Sections MII and MIII 
                %                                                     
                M(j,k) = M(j,k) + 1/2*sum(partials.gj.*partials.gk./(outg).^2);
              end
            end
          end
        end
      end
    end
    % based on M defined above, the convariance matrix is
    % defined as the inverse of M :
    % COV = M^-1
    s.ml.M = M;
    s.degrees_of_freedom     = num_measurements.all - length(parameters);
    s.covariance             = inv(M);
    s.num_observations       = num_measurements.all;
    s.objective              = objective;
    s.aic                    = 2*objective + 2*length(parameters);
    s.bic                    = 2*objective + log(s.num_observations)*length(parameters);
  end



  % calculating the cv% and confidence intervals
  for i=1:length(parameters)
      s.coefficient_of_variation(i) = 100.*sqrt(s.covariance(i,i))./parameters(i);
      s.confidence_interval.lower_bound(i) = parameters(i) - sqrt(s.covariance(i,i))*tinv(.975,s.degrees_of_freedom);
      s.confidence_interval.upper_bound(i) = parameters(i) + sqrt(s.covariance(i,i))*tinv(.975,s.degrees_of_freedom);
  end

