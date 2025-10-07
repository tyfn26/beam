% ---- Replace the previous element loop with this debug+auto-fix loop ----

for ie=1:nElems
    elemNodes = elems(ie,:);
    Xe = nodes(elemNodes, :); % 8x3

    % Quick duplicate-node check (zero-volume risk)
    dups = false;
    for p=1:8
        for q=p+1:8
            if norm(Xe(p,:) - Xe(q,:)) < 1e-12
                fprintf('Element %d has duplicate nodes: %d and %d (coords identical)\n', ie, p, q);
                dups = true;
            end
        end
    end
    if dups
        error('Element %d has duplicate nodes -> non-positive/zero Jacobian likely.', ie);
    end

    % FUNCTION: test detJ at 2x2x2 gauss for a given node ordering
    function [minDetJ, dets] = min_detJ_for_order(node_order)
        Xe_try = nodes(node_order, :);
        dets = zeros(8,1);
        k = 0;
        for ii=1:2
            for jj=1:2
                for kk=1:2
                    xi = gp(ii); eta = gp(jj); zeta = gp(kk);
                    dN_dxi = 1/8*[
                        [ -(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta) ];
                        [  (1-eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1-eta) ];
                        [  (1+eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1+eta) ];
                        [ -(1+eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1+eta) ];
                        [ -(1-eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1-eta) ];
                        [  (1-eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1-eta) ];
                        [  (1+eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1+eta) ];
                        [ -(1+eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1+eta) ];
                    ];
                    J = dN_dxi' * Xe_try;
                    dets(k+1) = det(J);
                    k = k+1;
                end
            end
        end
        minDetJ = min(dets);
    end

    % Check original ordering
    [minDet_orig, dets_orig] = min_detJ_for_order(elemNodes);
    if minDet_orig > 1e-12
        % good: proceed normally with original ordering
        chosenOrder = elemNodes;
        usedFix = 'none';
    else
        % try a set of reasonable permutations/fixes
        triedOrders = {};
        fixes = {};
        % 1) reverse top face (swap order of 5..8)
        ord1 = elemNodes;
        ord1(5:8) = ord1([8 7 6 5]);
        triedOrders{end+1} = ord1; fixes{end+1} = 'reverse top face (5..8)';
        % 2) reverse bottom face (swap 1..4)
        ord2 = elemNodes;
        ord2(1:4) = ord2([4 3 2 1]);
        triedOrders{end+1} = ord2; fixes{end+1} = 'reverse bottom face (1..4)';
        % 3) swap top and bottom faces (1<->5,2<->6,3<->7,4<->8)
        ord3 = elemNodes;
        ord3([1 2 3 4 5 6 7 8]) = elemNodes([5 6 7 8 1 2 3 4]);
        triedOrders{end+1} = ord3; fixes{end+1} = 'swap top/bottom faces';
        % 4) reverse node winding (reverse all)
        ord4 = elemNodes([2 1 4 3 6 5 8 7]); % simple swap pairs to change winding
        triedOrders{end+1} = ord4; fixes{end+1} = 'reverse winding pairs';
        % 5) full reversed order
        ord5 = elemNodes([8 7 6 5 4 3 2 1]);
        triedOrders{end+1} = ord5; fixes{end+1} = 'full reverse';
        % 6) última fallback: swap nodes 2 and 4 (common fix)
        ord6 = elemNodes;
        tmp = ord6(2); ord6(2) = ord6(4); ord6(4) = tmp;
        triedOrders{end+1} = ord6; fixes{end+1} = 'swap 2 & 4';

        chosenOrder = [];
        usedFix = 'none';
        for t=1:length(triedOrders)
            ordTry = triedOrders{t};
            [minDet_try, dets_try] = min_detJ_for_order(ordTry);
            if minDet_try > 1e-12
                chosenOrder = ordTry;
                usedFix = fixes{t};
                break;
            end
        end

        if isempty(chosenOrder)
            % None of the simple fixes worked -> print diagnostics and error
            fprintf('Element %d: all simple reorder attempts failed. minDet(original)=%.3e\n', ie, minDet_orig);
            fprintf('Node coordinates (original ordering):\n');
            for a=1:8
                fprintf(' %d : (%.6g, %.6g, %.6g)\n', a, Xe(a,1), Xe(a,2), Xe(a,3));
            end
            error('Cannot fix element %d automatically. Inspect mesh/element connectivity.', ie);
        else
            fprintf('Element %d: fixed orientation using: %s (minDet after fix > 0)\n', ie, usedFix);
            % reflect the change in global elems array so future assembly uses fixed order
            elems(ie,:) = chosenOrder;
            Xe = nodes(chosenOrder, :);
        end
    end

    % Now with chosenOrder (Xe) compute Ke and Me as before
    % (use the previously shown B, NN assembly code but with nodes= chosenOrder)
    Ke = zeros(24,24);
    Me = zeros(24,24);
    % 2x2x2 Gauss integration
    for i=1:2
        for j=1:2
            for k=1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                w = gw(i)*gw(j)*gw(k);
                N = 1/8*[ (1-xi)*(1-eta)*(1-zeta);
                         (1+xi)*(1-eta)*(1-zeta);
                         (1+xi)*(1+eta)*(1-zeta);
                         (1-xi)*(1+eta)*(1-zeta);
                         (1-xi)*(1-eta)*(1+zeta);
                         (1+xi)*(1-eta)*(1+zeta);
                         (1+xi)*(1+eta)*(1+zeta);
                         (1-xi)*(1+eta)*(1+zeta)];
                dN_dxi = 1/8*[
                    [ -(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta) ];
                    [  (1-eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1-eta) ];
                    [  (1+eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1+eta) ];
                    [ -(1+eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1+eta) ];
                    [ -(1-eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1-eta) ];
                    [  (1-eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1-eta) ];
                    [  (1+eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1+eta) ];
                    [ -(1+eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1+eta) ];
                ];
                J = dN_dxi' * Xe;
                detJ = det(J);
                if detJ <= 0
                    error('Element %d still has non-positive detJ after attempted fix (detJ=%.3e).', ie, detJ);
                end
                invJ = inv(J);
                dN_dX = (invJ * dN_dxi')';
                B = zeros(6,24);
                for a=1:8
                    idx = (a-1)*3 + (1:3);
                    dNx = dN_dX(a,1); dNy = dN_dX(a,2); dNz = dN_dX(a,3);
                    B(:, idx) = [ dNx,    0,    0;
                                   0,   dNy,    0;
                                   0,    0,   dNz;
                                   dNy,  dNx,   0;
                                   0,    dNz,  dNy;
                                   dNz,  0,    dNx ]';
                end
                Ke = Ke + (B' * D * B) * detJ * w;
                NN = zeros(3,24);
                for a=1:8
                    Ni = N(a);
                    NN(:, (a-1)*3 + (1:3)) = Ni * eye(3);
                end
                Me = Me + rho * (NN' * NN) * detJ * w;
            end
        end
    end

    % Assemble into global
    gdof = zeros(24,1);
    for a=1:8
        gn = elems(ie,a); % note: elems updated if fix applied
        gdof((a-1)*3 + (1:3)) = (gn-1)*3 + (1:3);
    end
    K(gdof, gdof) = K(gdof, gdof) + Ke;
    M(gdof, gdof) = M(gdof, gdof) + Me;
end
% -------------------------------------------------------------------------
