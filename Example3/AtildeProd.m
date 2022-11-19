function prodA = AtildeProd(Atilde,baslangic,bitis)
    prodA=eye(size(Atilde{baslangic}));
    for i=baslangic:1:bitis
        prodA=Atilde{i}*prodA;
    end
end

