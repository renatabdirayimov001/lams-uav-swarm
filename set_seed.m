function set_seed(seed)
%SET_SEED  Seed the random generators in both MATLAB and GNU Octave.
if exist('OCTAVE_VERSION','builtin')
    rand('state', seed); randn('state', seed); %#ok<RAND>
else
    rng(seed, 'twister');
end
end
