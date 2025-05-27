
function out = option_input(in,change)

for k = 1:size(in,1)                 % C
    for n = 1:size(in,2)             % multiple
        for m = 1:size(in{k,n}.Z,2)  % soc
            p = in{k,n}.p;
switch change
    case 'cs_bar'
        out{k,n}{m} = [     in{k,n}.vari{m}.x                     in{k,n}.vari{m}.stoich];
    case 'cs_bar_neg'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(1:p.nn,:)           in{k,n}.vari{m}.stoich(1:p.nn,:)];
    case 'cs_bar_pos'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(p.nn+p.ns+1:end,:)  in{k,n}.vari{m}.stoich(p.nn+p.ns+1:end,:)];
    case 'ce'
        out{k,n}{m} = [     in{k,n}.vari{m}.x                     in{k,n}.vari{m}.ce_dim];
    case 'phis'
        out{k,n}{m} = [     in{k,n}.vari{m}.x                     in{k,n}.vari{m}.phis];
    case 'phis_neg'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(1:p.nn,:)           in{k,n}.vari{m}.phis(1:p.nn,:)];
    case 'phis_pos'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(p.nn+p.ns+1:end,:)  in{k,n}.vari{m}.phis(p.nn+p.ns+1:end,:)];
    case 'phie'
        out{k,n}{m} = [     in{k,n}.vari{m}.x                     in{k,n}.vari{m}.phie];
    case 'j'
        out{k,n}{m} = [     in{k,n}.vari{m}.x                     in{k,n}.vari{m}.j_dim];
    case 'j_neg'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(1:p.nn,:)           in{k,n}.vari{m}.j_dim(1:p.nn,:)];
    case 'j_pos'
        out{k,n}{m} = [     in{k,n}.vari{m}.x(p.nn+p.ns+1:end,:)  in{k,n}.vari{m}.j_dim(p.nn+p.ns+1:end,:)];
    case 'Z_neg'
        out{k,n}{m} = [real(in{k,n}.Z{m}.DFN_neg)                 imag(in{k,n}.Z{m}.DFN_neg)];
    case 'Z_pos'
        out{k,n}{m} = [real(in{k,n}.Z{m}.DFN_pos)                 imag(in{k,n}.Z{m}.DFN_pos)];
    case 'Z_sep'
        out{k,n}{m} = [real(in{k,n}.Z{m}.DFN_sep)                 imag(in{k,n}.Z{m}.DFN_sep)];
    case 'Z_cell'
        out{k,n}{m} = [real(in{k,n}.Z{m}.DFN_cell)                imag(in{k,n}.Z{m}.DFN_cell)];
    case 'ZDs_neg'
        out{k,n}{m} = [real(in{k,n}.Z{m}.Ds_neg)                  imag(in{k,n}.Z{m}.Ds_neg)];
    case 'ZDs_pos'
        out{k,n}{m} = [real(in{k,n}.Z{m}.Ds_pos)                  imag(in{k,n}.Z{m}.Ds_pos)];
    case 'ZDe_neg'
        out{k,n}{m} = [real(in{k,n}.Z{m}.De_neg)                  imag(in{k,n}.Z{m}.De_neg)];
    case 'ZDe_pos'
        out{k,n}{m} = [real(in{k,n}.Z{m}.De_pos)                  imag(in{k,n}.Z{m}.De_pos)];
end
        end
    end
end

end
