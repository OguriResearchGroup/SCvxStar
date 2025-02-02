function obj = load_scp(filename)
    % LOADSCP Load a SCP problem from a .mat file
    S = load(filename, 'obj');
    obj = S.obj;
end

