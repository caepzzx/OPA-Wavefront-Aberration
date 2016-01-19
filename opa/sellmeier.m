function n=sellmeier(lambda,material,type,var)
%SELLMEIER 此处显示有关此函数的摘要
%   此处显示详细说明
lambda=lambda*1e6;
switch material
    case 'BBO2'
%         switch type
%             case 'o'
%                 n=sqrt(2.7359+.01878./(lamda.^2-.01822)-.01354*lamda.^2+6.081e-4*lamda.^4-6.740e-5*lamda.^6);
%             case 'e'
%                 n=sqrt(2.3753+.01224./(lamda.^2-.01667)-.01516*lamda.^2+5.716e-4*lamda.^4-6.305e-5*lamda.^6);
%             case 't'
%                 no=sqrt(2.7359+.01878./(lamda.^2-.01822)-.01471*lamda.^2+6.081e-4*lamda.^4-6.740e-5*lamda.^6);
%                 ne=sqrt(2.3753+.01224./(lamda.^2-.01667)-.01516*lamda.^2+5.716e-4*lamda.^4-6.305e-5*lamda.^6);
%                 n=1./sqrt(sind(var).^2./ne.^2+cosd(var).^2./no.^2);
%         end
                switch type
            case 'o'
        A=2.7359;
        B=0.01878;
        C=0.01822;
        D=0.01471;
        E=0.0006081;
        F=0.0000674;
        
%文献<broad bandwidth parametric amplification in the visible femtosecond ...
...experiments and simulations> G.M.Gale 用的是下面的参数        
%         A=2.7405;
%         B=0.0184;
%         C=0.0179;
%         D=0.0155;
        
        n=sqrt(A+B./(lambda.^2-C)-D.*lambda.^2+E.*lambda.^4-F*lambda.^6);
%         n=sqrt(A+B./(lambda.^2-C)-D.*lambda.^2);
            case 'e'
        A=2.3753;
        B=0.01224;
        C=0.01667;
        D=0.01627;
        E=0.0005716;
        F=0.00006305;

 %文献<broad bandwidth parametric amplification in the visible femtosecond ...
...experiments and simulations> G.M.Gale 用的是下面的参数     
%         A=2.3730;
%         B=0.0128;
%         C=0.0156;
%         D=0.0044;
        
       n=sqrt(A+B./(lambda.^2-C)-D.*lambda.^2+E.*lambda.^4-F*lambda.^6);
%        n=sqrt(A+B./(lambda.^2-C)-D.*lambda.^2);
%BBO 系数取自《非线性光学晶体  一份完整的总结》
                end
        
    case 'BBO'
        
        switch type
            case 'o'
                n=sqrt(2.7359+0.01878./(lambda.^2-0.01822)-0.01354.*lambda.^2);
            case 'e'
                n=sqrt(2.3753+0.01224./(lambda.^2-0.01667)-0.01516.*lambda.^2);
            case 't'%相位匹配角为theta,对应的非寻常光的折射率
                no=sqrt(2.7359+0.01878./(lambda.^2-0.01822)-0.01354.*lambda.^2);
                ne=sqrt(2.3753+0.01224./(lambda.^2-0.01667)-0.01516.*lambda.^2);
                n=1./sqrt(sind(var).^2./ne.^2+cosd(var).^2./no.^2);


        end
        
    case 'LBO'
        nx=sqrt(2.4542+.01125./(lambda.^2-.01135)-.01388*lambda.^2);
        ny=sqrt(2.5390+.01277./(lambda.^2-.01189)-.01849*lambda.^2+4.3025e-5*lambda.^4-2.9131e-5*lambda.^6);
        nz=sqrt(2.4542+.01310./(lambda.^2-.01223)-.01862*lambda.^2+4.5778e-5*lambda.^4-3.2526e-5*lambda.^6);
   
    case 'PPLN'
        switch type
            case 'e'
% 1
%         f=(var-24.5).*(var+570.82);
%         a1=5.319725;
%         a2=.09147285;
%         a3=.3165008;
%         a4=100.2028;
%         a5=11.37639;
%         a6=1.497046e-2;
%         b1=4.753469e-7;
%         b2=3.310965e-8;
%         b3=2.760513e-5;
%         n=sqrt(a1+b1.*f+(a2+b2*f)./(lambda.^2-a3.^2)+(a4+b3.*f)./(lambda.^2-a5.^2)-a6.*lambda.^2);

% 2 中红外ppln光参量振荡技术研究 魏星斌 page13
%                 a1=5.756;
%                 a2=0.0983;
%                 a3=0.2020;
%                 a4=189.32;
%                 a5=12.52;
%                 a6=1.32e-2;
%                 b1=2.860e-6;
%                 b2=4.700e-8;
%                 b3=6.113e-8;
%                 b4=1.156e-4;
%                 f=(var-24.5).*(var+570.82);
%                 n=sqrt(a1+b1.*f+(a2+b2.*f)./(lambda.^2-(a3+b3.*f).^2)+(a4+b4.*f)./(lambda.^2-a5.^2)-a6.*lambda.^2);


%3 <<Sellmeier and thermo -optic dispersion formula s for the extra ordinary ray of 5 mol.% MgO-doped congruent LiNbO3 in the visible,infrared,and terahertz regions>>Nobuhiro Umemura
% n_temp=sqrt(4.54514+.096471./(lambda.^2-.043763)-.021502*lambda.^2).*(lambda<=2.7)+sqrt(24.6746+.0456./(lambda.^2-2.280)+19166.65./(lambda.^2-953.52)+1.0103./(lambda.^2-45.86)).*(lambda>2.7);
% n=n_temp+1e-5*((var-20)+.00138*(var-20).^2).*(.4175./lambda.^3-.6643./lambda.^2+.9036./lambda+3.5332-.0744.*lambda);

%4   0.4-5um 20-400°C 宽带准相位匹配飞秒光参量放大的研究 张海清 page25
                a1=5.35583;
                a2=.100473;
                a3=.20692;
                a4=100.0;
                a5=11.34927;
                a6=1.5334e-2;
                b1=4.629e-7;
                b2=3.862e-8;
                b3=-.89e-8;
                b4=2.657e-5;
                f=(var-24.5).*(var+570.5);
                n=sqrt(a1+b1.*f+(a2+b2.*f)./(lambda.^2-(a3+b3.*f).^2)+(a4+b4.*f)./(lambda.^2-a5.^2)-a6.*lambda.^2);
                








            case 'o'
                a1=5.653;
                a2=.1185;
                a3=.2091;
                a4=89.61;
                a5=10.85;
                a6=1.97e-2;
                b1=7.941e-7;
                b2=3.134e-8;
                b3=-4.641e-9;
                b4=-2.188e-6;
                f=(var-24.5).*(var+570.82);
                n=sqrt(a1+b1.*f+(a2+b2.*f)./(lambda.^2-(a3+b3.*f).^2)+(a4+b4.*f)./(lambda.^2-a5.^2)-a6.*lambda.^2);
                
        end
%         switch type
%             case 'o'
%                 n=sqrt(1+2.4272*lamda.^2./(lamda.^2-.01478)+1.4617*lamda.^2./(lamda.^2-.05612)+9.6536*lamda.^2./(lamda.^2-371.216));
%             case 'e'
%                 n=sqrt(1+2.2454*lamda.^2./(lamda.^2-.01242)+1.3005*lamda.^2./(lamda.^2-.05313)+6.8972*lamda.^2./(lamda.^2-331.33));
%             case 't'
%                 no=sqrt(1+2.4272*lamda.^2./(lamda.^2-.01478)+1.4617*lamda.^2./(lamda.^2-.05612)+9.6536*lamda.^2./(lamda.^2-371.216));
%                 ne=sqrt(1+2.2454*lamda.^2./(lamda.^2-.01242)+1.3005*lamda.^2./(lamda.^2-.05313)+6.8972*lamda.^2./(lamda.^2-331.33));
%                 n=1./sqrt(sind(theta).^2./ne.^2+cosd(theta).^2./no.^2);
%         end
        
    case 'FS'
        n=sqrt(1+.6961663*lambda.^2./(lambda.^2-.0684043^2)+.4079426*lambda.^2./(lambda.^2-.1162414^2)+.8974794*lambda.^2./(lambda.^2-9.896161^2));
    
    case 'SILICON'
        n=sqrt(1+.6961663.*lambda.^2./(lambda.^2-.0684043^2)+.4079426*lambda.^2./(lambda.^2-.1162414^2)+.8974794*lambda.^2./(lambda.^2-9.896161^2));

    case 'BK7'
        n=sqrt(1+1.03961212*lambda.^2./(lambda.^2-.00600069867)+.231792344*lambda.^2./(lambda.^2-.0200179114)+1.01046945*lambda.^2./(lambda.^2-103.560653));

end

