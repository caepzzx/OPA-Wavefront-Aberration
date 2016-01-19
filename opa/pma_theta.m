function [theta]=pma_theta(N,alpha,lamdas,lamdap,material,type )
%type:匹配方式sip对应着三种匹配方式ｏｏｅ，ｏｅｅ，ｅｏｅ
%condition: 共线方式
%
%
c=3e8;
switch type
    case 'ooe'
        % os+oi->ep
        no_s=sellmeier(lamdas,material,'o');
        lamdai=1./(1./lamdap-1./lamdas);
        no_i=sellmeier(lamdai,material,'o');
        no_p=sellmeier(lamdap,material,'o');
        ne_p=sellmeier(lamdap,material,'e');
               
        ks=2*pi*no_s./lamdas;
        ki=2*pi*no_i./lamdai;
        gamma=asind(ks*sind(alpha)./ki);%信号光与闲置光动量垂直分量相等
        kp=ks.*cosd(alpha)+ki.*cosd(gamma);%动量水平分量守恒
        theta=asind(ne_p./kp.*sqrt(((2*pi./lamdap).^2.*no_p.^2-kp.^2)./(no_p.^2-ne_p.^2)));%参考文献current applied physics 12(2012)648-653 page649 公式(5)


    case 'oee'
        %os+ei->ep
        the=linspace(0,90,1000);
        no_s=sellmeier(lamdas,material,'o');
        ws=2*pi*c./lamdas;
        lamdai=1./(1./lamdap-1./lamdas);
        no_i=sellmeier(lamdai,material,'o');
        ne_i=sellmeier(lamdai,material,'e');
        wi=2*pi*c./lamdai;
        no_p=sellmeier(lamdap,material,'o');
        ne_p=sellmeier(lamdap,material,'e');
        wp=2*pi*c/lamdap;
        theta=zeros(1,N); 
        
        for i=1:N
            ne_it=1./sqrt(sind(the).^2./ne_i(i)^2+cosd(the).^2./no_i(i).^2);
            ne_pt=1./sqrt(sind(the).^2./ne_p^2+cosd(the).^2./no_p.^2);
            [q,p]=min(abs(no_s(i)*ws(i)+ne_it*wi(i)-ne_pt*wp));
            theta(i)=the(p);
        end
   

     

    case 'eoe'
        % es+oi->ep
        the=linspace(0,90,1000);
        no_s=sellmeier(lamdas,material,'o');
        ne_s=sellmeier(lamdas,material,'e');
        ws=2*pi*c./lamdas;
        lamdai=1./(1./lamdap-1./lamdas);
        no_i=sellmeier(lamdai,material,'o');
        wi=2*pi*c./lamdai;
        no_p=sellmeier(lamdap,material,'o');
        ne_p=sellmeier(lamdap,material,'e');
        wp=2*pi*c/lamdap;
        theta=zeros(1,N);

        for i=1:100
            ne_st=1./sqrt(sind(the).^2./ne_s(i)^2+cosd(the).^2./no_s(i).^2);
            ne_pt=1./sqrt(sind(the).^2./ne_p^2+cosd(the).^2./no_p.^2);
            [q,p]=min(abs(ne_st*ws(i)+no_i(i)*wi(i)-ne_pt*wp));
            theta(i)=the(p);
        end

end

