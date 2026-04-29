
for i=1:2*s+1
    for j=1:2*s+1
    
        x(i,j)=P_trial(1,(j-1)*(2*s+1)+i);
        y(i,j)=P_trial(2,(j-1)*(2*s+1)+i);
    end
end
for i=1:s+1
    for j=1:s+1
    
        xp(i,j)=P(1,(j-1)*(s+1)+i);
        yp(i,j)=P(2,(j-1)*(s+1)+i);
    end
end

u1 = ul(1:length(P_trial));          % 速度x分量（示例：正弦分布）
u2 = ul(length(P_trial)+1:2*length(P_trial));          % 速度y分量（示例：余弦分布）
p = ul(2*length(P_trial)+1:end);  % 压力（示例：中心高、边缘低的高斯分布）


plot_flow_component(x, y, xp, yp, u1, u2, p,  'channel');
