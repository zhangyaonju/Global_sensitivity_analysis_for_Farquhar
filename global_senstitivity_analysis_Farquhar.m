%% global sensitivity analysis
% the variance based GSA.
% Sobol's technique
% for the Farquhar model coupled with a big leaf model

% range of para
pmin = [0, 50, 30000, 0.2, 10000, 10, 10000, 10, 30000, 0.7, 1, 2, 0.2];
pmax = [0.5, 600, 150000, 0.5, 60000, 300, 100000, 150, 100000, 0.9, 5, 30, 10];
npara = length(pmin);

N = 5000;
k = npara;
onemat = fnc_getSobolSequence(k, N);
P = tomaxmin(pmin,pmax,onemat);
onemat2 = fnc_getSobolSequence(k+3, N+5000);
Q = tomaxmin(pmin,pmax,onemat2(5001:(N+5000),4:k+3));
%%plot(P(:,3),Q(:,3),'.');

%% get the site information
%sysdir = '/data/ifs/users/yzhang/PROJECT/MCMC/';
sysdir = 'V:/users/yzhang/PROJECT/MCMC/';
text = importdata(strcat(sysdir,'site_info.txt'),'');
site_id = char(text(2:end));
%% read data from each site-year
listing = dir(strcat(sysdir,'input_data_with_uncert_site_year/*.csv'));

siteid = [172,244,254,278];

for o = 1:4  
    sid = siteid(o);
    file = listing(sid).name;
    mkdir(strcat(sysdir,'sensitivity_test/',file(1:6)));
    input = csvread(strcat(sysdir,'input_data_with_uncert_site_year/',file),1,1);
    
    I = input(:,2)*2.148;
    Tk = input(:,3)+273.15;
    D = input(:,4);
    Ca = input(:,6);
    LAI = input(:,7);
    %GPP = input(:,8);
    %GPPuncert = input(:,9);
    GPPsimuP = zeros(N,1);
    GPPsimuQ = zeros(N,1);
    for n = 1:N
        GPPsimuP(n) = nansum(Farquhar_model(P(n,:),I,Tk,D,Ca,LAI));
        GPPsimuQ(n) = nansum(Farquhar_model(Q(n,:),I,Tk,D,Ca,LAI));
    end
    y=Farquhar_model(P(882,:),I,Tk,D,Ca,LAI);
    %2 digit month as char
    GPPsimurepP = zeros(N,k);
    GPPsimurepQ = zeros(N,k);
    for para = 1:k
        for n = 1:N
            GPPsimurepP(n,para) = nansum(Farquhar_model(replacepar(P(n,:),Q(n,:),para),I,Tk,D,Ca,LAI));
            GPPsimurepQ(n,para) = nansum(Farquhar_model(replacepar(Q(n,:),P(n,:),para),I,Tk,D,Ca,LAI));
        end
        %Vhat(para) = 1/(N-1)* nansum(GPPsimuP.*GPPsimurepP(:,para))-fnotsqu;
    end;
    
    %% calculate the stats
    S = zeros(4991,13);    
    for m = 10:5000
        fnotsqu = mean(GPPsimuP(1:m).*GPPsimuQ(1:m));
        varYnot = 1/2/(m-1)*(nansum(power(GPPsimuP(1:m),2))+nansum(power(GPPsimuQ(1:m),2)))-fnotsqu;
        for i = 1:k
            Vhat(i) = 1/(m-1)* nansum(GPPsimuP(1:m).*(GPPsimurepQ(1:m,i)))-fnotsqu;
        end
        
        for i = 1:k
            for j = (i+1):k
                Vhatij(i,j) = 1/(m-1)*nansum(GPPsimurepQ(1:m,i).*GPPsimurepP(1:m,j))-fnotsqu-Vhat(i)-Vhat(j);
            end
            Vhatnoti(i) = 1/(m-1)*nansum(GPPsimuQ(1:m).*GPPsimurepQ(1:m,i))-fnotsqu;
        end
        
        S(m-9,1:13) = Vhat/varYnot;
    end  
    
    semilogy(S)
    
    filename = strcat(sysdir,'sensitivity_test/',file(1:6),'/',...
        file(1:11),'.csv');
    csvwrite(filename,S);
   
end


