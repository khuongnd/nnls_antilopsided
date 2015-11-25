function his =  FCD(W,x,opt)

Q = full(W'*W);
q = full(-x'*W)';

startTime = tic;
[~,r]=size(W);
[m,n]=size(x);

alpha=zeros(r,n);
w=zeros(m,n);

maxIter=opt.maxIter;

lambda=0;

for k=1:r,
    normw(k)=norm(W(:,k))^2;
end

Wtx=W'*x;

totalMinusTime = 0;
minusStart = tic;
Grad = Q*alpha+q;
his.x = alpha;
BB = (alpha > 0 | Grad < 0);
        
his.obj = norm(x - W*alpha, 'fro').^2;
his.d_barf = sum(Grad(BB).^2);
his.time = toc(startTime);
his.x = alpha;
totalMinusTime = totalMinusTime + toc(minusStart);

for i=1:maxIter,
    
    for j=randperm(r),
        G = W(:,j)'*w - Wtx(j,:)  + lambda ;
        ind=find(alpha(j,:)==0);
        pg = G;
        pg(ind)=min(G(ind),0);
  
        ind=find(abs(pg)~=0);
        if prod(size(ind))~=0
            alphatmp=alpha(j,ind);
            alpha(j,ind) = max(alpha(j,ind) - G(1,ind)./(normw(j)),0);
            %w(:,ind) = w(:,ind) + W(:,j)*ones(1,size(ind,2))*sparse(diag(((alpha(j,ind) - alphatmp))));
            w(:,ind) = w(:,ind) + W(:,j)*(alpha(j,ind) - alphatmp);
        end  
    end
    
    minusStart = tic;
    Grad = Q*alpha+q;
    his.x = alpha;
    BB = (alpha > 0 | Grad < 0);
        
    his.obj = [his.obj norm(x - W*alpha, 'fro').^2];
    his.d_barf = [his.d_barf sum(Grad(BB).^2)];
    totalMinusTime = totalMinusTime + toc(minusStart);
    his.time = [his.time (toc(startTime) - totalMinusTime)];
end  
his.finalObj = norm(x - W*alpha, 'fro')^2/2.0;
%h=alpha;
%time=toc;
end

