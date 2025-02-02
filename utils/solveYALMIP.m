function [sol, solve_flag] = solveYALMIP(constraints, objective, yalmip_options)
    sol = optimize(constraints, objective, yalmip_options);
    solve_flag = sol.problem;
    switch sol.problem
        case 0 % Solved
            return
        case 1 % Infeasible. resolve problem without objective 
            fprintf('Subproblem not solved: %s\n', sol.info)
            feasibility_problem = optimize(constraints, [], yalmip_options);
            if feasibility_problem.problem == 0
                fprintf('Feasibility problem solved. There may be a problem in the objective. \n')                        
            else
                fprintf('Feasibility problem not solved. The subproblem is infeasible. \n')
            end
            return
        case 4 % Numerical problems
            fprintf('Numerical problems: MOSEK Error Code %i\n', sol.solveroutput.r)
        otherwise
            fprintf('Subproblem not solved: %s\n', sol.info)
    end

end