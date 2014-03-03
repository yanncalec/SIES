classdef Frame
% Class of general frame operator. 
% 
% A (wavelet) frame is constituted by a shift-invariant space called approximation space (generated
% by translating a mother scaling function), and several shift-invariant spaces called detail spaces
% (generated by translating, dilating (and rotating) a mother wavelet function). Both of them are
% generated using a single mother function.
%     
% Convention: in this class and other related functions, the scale index ranging from 1 to nbScl is
% denoted by j, while the true scale value ranging from Jmin to Jmax is denoted by the capital J. We
% call j the scale index and J the scale number (or simply the scale). Similarly, the position
% number n (the integer n of translation in 2^J*n, n\in\Z^2) differs from the position index ranging
% from 1 to maximum size of an array.
%
% For ApprxSpace and DetailSpace: the small index corresponds to the coarse scale, e.g.,
% ApprxSpace{1} is the coarsest approximation space, and ApprxSpace{end} is the finest one.
    
    properties(SetAccess=protected)
        Rsl % spatial resolution of the imaging system
        
        Jmin % the finest scale of detailed spaces
        Jmax % the coarsest scale of detailed spaces
        nbScl % number of detail spaces = Jmax-Jmin+1
        
        sclFct=2 % sclFct>1 is the scaling factor between detail spaces
        splStep=1 % sampling step for translation in detail spaces
        apprx_splStep=1 % sampling step in the approximation space
        
        epsilon=0 % truncation threshold for wavelet's support

        ApprxSpace % shift-invariant spaces of approximation, instance of wavelet.ShiftInvSpace class
        DetailSpace % shift-invariant spaces of details, a cell, instance of wavelet.ShiftInvSpace class
        
        ROI % square ROI, instance of wavelet.ROI class

        nbWvl % total number of functions in the frame, sum of nbDetail and nbApprx
        nbDetail % total number of wavelets (detail spaces) in all scales, directions and positions
        nbApprx % total number of scaling functions (apprx space) in all positions        
        
        wvlname % the name of the wavelet frame
    end

    methods                
        function val = get.nbScl(obj)
            val = obj.Jmax-obj.Jmin+1;
        end
        
        function val = get.nbWvl(obj)
            val = obj.nbApprx + obj.nbDetail;
        end
        
        function val = get.nbApprx(obj)
            val = prod(obj.ApprxSpace{1}.dim);
        end
        
        function val = get.nbDetail(obj)
            val = 0;

            for j=1:obj.nbScl
                for k=1:length(obj.DetailSpace(j,:))
                    val = val + prod(obj.DetailSpace{j,k}.dim);
                end
            end
        end        
        
        function val = nbWvl_lowscls(obj, nbScl)
        % Number of wavelets in the first nbScl scales (approximation space included)
            val = obj.nbApprx;
            for j=1:nbScl
                val = val + prod(obj.DetailSpace{j,1}.dim) + prod(obj.DetailSpace{j,2}.dim) + prod(obj.DetailSpace{j,3}.dim);
            end
        end
    
    end
    
    methods(Abstract)        
        X = analysis(obj, F) % analysis operator
        X = synthesis(obj, C) % synthesis operator        
        X = dual_analysis(obj, F) % analysis operator by dual frame
        X = dual_synthesis(obj, C) % synthesis operator by dual frame

        val = eval_apprx(obj, X) % evaluate the mother scaling function
        val = eval_detail(obj, X, J, k) % evaluate the mother wavelet function

        W = WPT_std(obj, D, freq) % Wavelet Polarization Tensor, standard form
        W = WPT_nonstd(obj, D, freq) % Wavelet Polarization Tensor, non standard form
        
        W = WPT_vec2wvl(obj, X) % convert a vector to a WPT structure
        X = WPT_wvl2vec(obj, W) % convert a WPT structure to a vector

        % Calculate the representation of 2D Green function of Laplacian with the frame, such
        % that the resynthesized function is close to the original Green function on
        % the ROI.
        val = decomp_Green2D(obj, X0)        
    end
    
    methods(Abstract, Static)
        val = active_position(Z0, Lx, Ly, J, k) % get the position of wavelets/scaling functions
                                                % contributing (intersection of their supports)
                                                % to the ROI
    end
    
end
