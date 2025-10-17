classdef EasyLM
    properties (SetAccess = private)
        Params_
        OptFunc_
        UpdateFunc_
        Fix_id_
        Unfix_id_
    end
    methods(Access = public)
        function obj = EasyLM(opt_func, update_func, varargin)
            obj.Params_ = inputParser;
            obj.Params_.KeepUnmatched = true;
            obj.Params_.addParameter('MaxIteration', 100);
            obj.Params_.addParameter('miu', 0.01);
            obj.Params_.addParameter('tolX', 1e-8);
            obj.Params_.addParameter('tolFun', 1e-8);
            obj.Params_.addParameter('tolOpt', 1e-10);
            obj.Params_.addParameter('SpecifyObjectiveGradient', false);
            obj.Params_.addParameter('CheckGradient', false);
            obj.Params_.addParameter('ParallelNumericalDiff', false);
            obj.Params_.addParameter('FixParameter', []);
            obj.Params_.addParameter('Debug', 0);
            obj.Params_.parse(varargin{:});
            obj.OptFunc_ = opt_func;
            obj.UpdateFunc_ = update_func;


        end


        function x_opt = solve(obj, x)

            if ~isempty(obj.Params_.Results.FixParameter)
                obj.Fix_id_ = obj.Params_.Results.FixParameter;
                obj.Unfix_id_ = setdiff((1:size(x, 2)), obj.Fix_id_);
            end


            miu = obj.Params_.Results.miu;
            tolX = obj.Params_.Results.tolX;
            tolFun = obj.Params_.Results.tolFun;
            tolOpt = obj.Params_.Results.tolOpt;


            iter = 0;
            nu = 2;
            sqrtEps=sqrt(eps);


            [F,J] = obj.Calculate_F_J(x);
            sqSumF = dot(F,F);
            Residual = sqSumF;
            if ~isempty(obj.Fix_id_)
                J = J(:,obj.Unfix_id_);
            end
            H = J'*J;
            JtF = J'*F;
            %obj.disp_local([num2str(iter) ' : ' num2str(sqSumF)]);


            while iter < obj.Params_.Results.MaxIteration
                %sum(x)
                iter = iter+1;

                p_ori = x(:,obj.Unfix_id_);
                % --- calculate update ---
                H_LM = H + sparse(1:length(p_ori),1:length(p_ori),miu,length(p_ori),length(p_ori));
                dp = -H_LM \ JtF;

                % --- calculate update ---
                dx = zeros(size(x));
                dx(:,obj.Unfix_id_) = dp;
                x_LM = obj.UpdateFunc_(x, dx);

                % --- calculate new FJ ---
                [F_LM,J_LM] = obj.Calculate_F_J(x_LM);


                if ~isempty(obj.Fix_id_)
                    J_LM = J_LM(:,obj.Unfix_id_);
                end
                sqSumF_LM = dot(F_LM,F_LM);
                JtF_LM = J_LM'*F_LM;


                if (norm(dp) < tolX*(sqrtEps + norm(x)))
                    obj.disp_local('Finished (tolX)');
                    break;
                end

                if (norm(JtF_LM,Inf) < tolOpt)
                    obj.disp_local('Finished (tolOpt)');
                    break;
                end

                if (abs(sqSumF_LM - sqSumF) <= tolFun*(sqrtEps + sqSumF) )
                    obj.disp_local('Finished (tolFun)');
                    break;
                end

                %varrho = (sqSumF_LM - sqSumF) / (dp' * JtF_LM);
                varrho = -(sqSumF_LM - sqSumF);

                % if (varrho > 0)
                %     x = x_LM;
                %     sqSumF = sqSumF_LM;
                %     JtF = JtF_LM;
                %     H = J_LM'*J_LM;
                % 
                % 
                %     obj.disp_local([num2str(iter) ' : ' num2str(sqSumF)]);
                %     Residual = [Residual, sqSumF];
                % 
                %     miu = miu * max(0.3333, 1-(2*varrho-1)^3);
                %     nu = 2;
                % else
                %     miu = miu * nu;
                %     nu = nu * 2;
                % end

                if (varrho > 0)
                    miu=miu*0.1;
                    x=x_LM;
                    sqSumF=sqSumF_LM;
                    %disp([num2str(iter) ' : ' num2str(sqSumF)]);
                    Residual = [Residual, sqSumF];
                    H = J_LM'*J_LM;
                    JtF = JtF_LM;
                else
                    miu=miu*10;
                end

            end

            if iter >= obj.Params_.Results.MaxIteration
                obj.disp_local('Finished (max_iter)');
            end

            x_opt=x_LM;
        end

    end
    methods(Access = private)


        function [F,J] = Calculate_F_J(obj, x)

            if obj.Params_.Results.SpecifyObjectiveGradient
                [F,J] = obj.OptFunc_(x);
                if obj.Params_.Results.CheckGradient
                    J_num = obj.NumericalDiff(x, F);
                    obj.disp_local(['Check Gradient: ', num2str(norm(J - J_num))]);
                end
            else
                [F] = obj.OptFunc_(x);
                J = obj.NumericalDiff(x, F);
            end

        end


        function disp_local(obj, message)
            if obj.Params_.Results.Debug
                display(message);
            end
        end



        function [J] = NumericalDiff(obj, x_opt, F)
            parameter_dim = size(x_opt, 2);
            J = zeros(size(F,1), parameter_dim);

            if obj.Params_.Results.ParallelNumericalDiff

            else
                for i = 1:parameter_dim
                    % --- forward ---
                    x_p = zeros(1, parameter_dim);
                    x_p(1, i) = 1e-6;
                    x_tmp = obj.UpdateFunc_(x_opt, x_p);
                    [F_p] = obj.OptFunc_(x_tmp);


                    % --- backward ---
                    x_n = zeros(1, parameter_dim);
                    x_n(1, i) = -1e-6;
                    x_tmp = obj.UpdateFunc_(x_opt, x_n);
                    [F_n] = obj.OptFunc_(x_tmp);

                    j_slice = (F_p - F_n) / 2e-6;
                    J(:, i) = j_slice';
                end
            end
        end

        function [J] = L_NumericalDiff(obj, x_opt, F)
            parameter_dim = size(x_opt, 2);
            J = zeros(size(F,1), parameter_dim);

            if obj.Params_.Results.ParallelNumericalDiff

            else
                for i = 1:parameter_dim
                    % --- forward ---
                    x_p = zeros(1, parameter_dim);
                    x_p(1, i) = 1e-6;
                    x_tmp = obj.UpdateFunc_(x_opt, x_p);
                    [F_p] = obj.OptFunc_(x_tmp);


                    % --- backward ---
                    x_n = zeros(1, parameter_dim);
                    x_n(1, i) = -1e-6;
                    x_tmp = obj.UpdateFunc_(x_opt, x_n);
                    [F_n] = obj.OptFunc_(x_tmp);

                    j_slice = (F_p - F_n) / 2e-6;
                    J(:, i) = j_slice';
                end
            end
        end



    end

end






