function [allPointMember,numNeiPoint,pointfamtot,ih,fail] = two_neighboornodes(coord,Geome)
delta = Geome.delta;
totnode = size(coord,1);
numNeiPoint = int32(zeros(totnode,1));
pointfamtot = int32(zeros(totnode+1, 1));
[Idx,~] = rangesearch(coord,coord,delta);
ih = int32(zeros(totnode*200,1));
for i = 1:totnode
    Idx{i,1} = sort(Idx{i,1}(2:end));
    %numfam(i,1) = size(Idx{i,1},2);
    numNeiPoint(i,1) = size(Idx{i,1},2);
    pointfamtot(i+1,1) = pointfamtot(i)+numNeiPoint(i,1);
    ih(pointfamtot(i,1)+1:pointfamtot(i+1,1),1) = i;
end
allPointMember = int32(([Idx{:}])');
ih = ih(ih>0);

TipL = Geome.TipL;
TipR = Geome.TipR;
XY1 = [TipL(:,1) TipL(:,2) TipR(:,1) TipR(:,2)];
%XY2 = [-60,10,-40,30;-20,-30,-20,-45;40,-40,60,0;-20,-20,20,20];
XY2 = [coord(ih,1),coord(ih,2),coord(allPointMember,1),coord(allPointMember,2)];
%line([TipL(:,1)';TipR(:,1)'],[TipL(:,2)';TipR(:,2)'],'Color','black','LineWidth',1)
%line([XY2(:,1)';XY2(:,3)'],[XY2(:,2)';XY2(:,4)'])
out = lineSegmentIntersect(XY2,XY1);
ax = sum(out.intAdjacencyMatrix,2);
index = ax > 0;
fail = ones(size(allPointMember,1),1);
fail(index) = 0;
%dmg = 1-accumarray(ih,fail.*pv(allPointMember,1))./accumarray(ih,pv(allPointMember,1));