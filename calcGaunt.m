function int = calcGaunt(l1,m1,l2,m2,LLoad0,C0,S0,w3js,Itab)
% int = calcGaunt(l1,m1,l2,m2,LLoad0,C0,S0,w3js,Itab)
%
% Calculates the integral of XlmXl'm' from the list of integrals of Xlm
% using the method of Gaunt 
%
% INPUT:
%
% l1, m1    degree, order of the first Xlm
% l2, m2    degree, order of the second Xlm
% LLoad0    Maximum degree of the loaded wigner0j database C0,S0
% C0,S0     Loaded wigner0j database
% w3js      Loaded wigner3j database
% Itab      List of integrals of the Xlm, from PAUL
%
% OUTPUT:
% 
% int       value of the calculated integral
%
% See also PAUL, KERNELBM, KERNELB
%
% Last modified by plattner-at-alumni.ethz.ch, 01/30/2012

if (l1<0 || l2<0)
    int=0;
else  
    % The ls in the Gaunt sum for
    ELL=(max(abs(l1-l2),abs(m1+m2)):(l1+l2))';     
    % Now calculate all Qs for the different ELLs              
    w0=zeroj(l1,l2,ELL,LLoad0,2,C0,S0);  
    [CC,oddperm,phasefix]=wignersort(ELL,l1,l2,-m1-m2,m1,m2);
    % Do the initial evaluation from the loaded variable
    % for threej
    wm=full(w3js(CC));
    % Fix the phase
    wm(oddperm)=wm(oddperm).*phasefix;
    % Now fix the triangle condition violations
    wm=wm.*triangle(repmat(l1,length(ELL),1),repmat(l2,length(ELL),1),ELL);
    % Now fix the order violations
    wm=wm.*~[l1<abs(m1) | l2<abs(m2) | ELL<abs(-m1-m2)];
    % Calculate the Gaunt coefficients
    Q=(-1)^(m1+m2)*(2*ELL+1).*wm .*w0';  
    % The Integral table Itab only contains the integrals for positive 
    % m1+m2. Therefore take the positive entries and multiply with 
    % (-1)^(m1+m2) if m1+m2 is negative (because Xl-m = (-1)^mXlm)
    indices=ELL.*(ELL+1)/2+abs(m1+m2)+1;
    sig=1;
    if(m1+m2 < 0)
       sig=(-1)^(m1+m2);
    end
    Itab=[Itab;zeros(2,size(Itab,2))];
    % Find the corresponding entries in Itab, multply them
    % with the Gaunt coefficients and sum   
    int=Q'*Itab(indices,:)*sig;
end
