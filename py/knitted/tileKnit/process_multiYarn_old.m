function process
    tiling = [5, 5];
    %tiling = [50, 50];

    outputObj = false;
    outputMts = false;
    plotResult = true;

    fin = fopen("slipstitchrib_2.3.txt", "rt");
    n = str2double(extractAfter(fgetl(fin), 11));
    fibers = cell(1, n);
    totPts = 0;
    for i = 1 : n
        m = fscanf(fin, "%d", 1);
        fibers{i} = zeros(3, m);
        for j = 1 : m
            fibers{i}(:, j) = fscanf(fin, "%f", 3);
        end
        totPts = totPts + m;
    end
    fclose(fin);
    totSeg = totPts - n;
    
    %{
    clf
    hold all
    for i = 1 : length(fibers)
        plot3(fibers{i}(1, :), fibers{i}(2, :), fibers{i}(3, :))    
        plot3(fibers{i}(1, :), fibers{i}(2, :) + 10, fibers{i}(3, :))
        plot3(fibers{i}(1, :) + 22.5, fibers{i}(2, :), fibers{i}(3, :)) 
    end
    hold off
    axis equal
    grid on
    xlabel('X')
    ylabel('Y')
    legend('Location', 'BestOutside')
    %}

    dx = 22.51; dy = 10;
    thres = 0.1;
    thres1 = 15.0;
    
    vtx0 = zeros(totPts, 3);
    seg0 = zeros(totSeg, 2);

    evtx0 = zeros(2*n, 3);
    eidx0 = zeros(2*n, 1);
    
    k = 0; t = 0;
    for i = 1 : n
        m = length(fibers{i});
        
        evtx0(2*i - 1, :) = fibers{i}(:, 1)';
        evtx0(2*i, :) = fibers{i}(:, end)';
        eidx0(2*i - 1) = k + 1;
        eidx0(2*i) = k + m;

        for j = 1 : m
            k = k + 1;
            vtx0(k, :) = [fibers{i}(1, j), fibers{i}(2, j), fibers{i}(3, j)];
            if j < m
                t = t + 1;
                seg0(t, :) = [k, k + 1];
            end
        end
    end
    
    totTile = tiling(1)*tiling(2);    
    dxMat = kron(ones(1, tiling(2)), (0 : tiling(1) - 1)');
    dyMat = kron(ones(tiling(1), 1), 0 : tiling(2) - 1);
    
    vtx = kron(ones(totTile, 1), vtx0) + ...
          kron(dxMat(:), kron(ones(totPts, 1), [dx, 0, 0])) + ...
          kron(dyMat(:), kron(ones(totPts, 1), [0, dy, 0]));
    seg = kron(ones(totTile, 1), seg0) + kron((0 : totTile - 1)', totPts*ones(totSeg, 2));
    
    evtx = kron(ones(totTile, 1), evtx0) + ...
           kron(dxMat(:), kron(ones(2*n, 1), [dx, 0, 0])) + ...
           kron(dyMat(:), kron(ones(2*n, 1), [0, dy, 0]));
    eidx = kron(ones(totTile, 1), eidx0) + kron((0 : totTile - 1)', totPts*ones(2*n, 1));
    totVtx = size(vtx, 1);
    totEVtx = size(evtx, 1);
  
    flag = zeros(1, totVtx);
    count = 0;
    ptCloud = pointCloud(evtx);
    for i = 1 : totEVtx
        i1 = eidx(i);
        if flag(i1) == 0
            [idx, dist] = findNearestNeighbors(ptCloud, evtx(i, :), 2);
            assert(idx(1) == i || idx(2) == i);
            if idx(1) == i
                k = 2;
            else
                k = 1;
            end
            if dist(k) < thres
                j = idx(k);
                j1 = eidx(j);
                assert(flag(j1) == 0)
                assert(norm(vtx(i1, :) - vtx(j1, :)) < thres)
                
                flag(j1) = i1;
                avg = 0.5*(vtx(i1, :) + vtx(j1, :));
                vtx(i1, :) = avg;
                vtx(j1, :) = avg;
                
                count = count + 1;
            end
        end
    end
    fprintf("%d vertex pairs merged.\n", count)
       
    gidx = zeros(1, totVtx);
    totVtx1 = 0;
    for i = 1 : totVtx
        if flag(i) == 0
            totVtx1 = totVtx1 + 1;
            gidx(i) = totVtx1;
        else
            gidx(i) = gidx(flag(i));
        end
    end
    
    vtx1 = zeros(totVtx1, 3);
    for i = 1 : totVtx
        vtx1(gidx(i), :) = vtx(i, :);
    end
    lower = [min(vtx1(:, 1)), min(vtx1(:, 2)), min(vtx1(:, 3))];
    upper = [max(vtx1(:, 1)), max(vtx1(:, 2)), max(vtx1(:, 3))];
    center = 0.5*(lower + upper);
    vtx1 = vtx1 - kron(ones(totVtx1, 1), center);   
    fprintf("Final AABB: [%.4f, %.4f, %.4f]--[%.4f, %.4f, %.4f]\n", ...
        min(vtx1(:, 1)), min(vtx1(:, 2)), min(vtx1(:, 3)), ...
        max(vtx1(:, 1)), max(vtx1(:, 2)), max(vtx1(:, 3)));
    
    lst = zeros(totVtx1, 3);
    for i = 1 : size(seg, 1)
        u = gidx(seg(i, 1));
        v = gidx(seg(i, 2));
        assert(u ~= v && u > 0 && v > 0);
        
        assert(lst(u, 1) < 2 && lst(v, 1) < 2)
        lst(u, 1) = lst(u, 1) + 1;
        lst(u, lst(u, 1) + 1) = v;
        lst(v, 1) = lst(v, 1) + 1;
        lst(v, lst(v, 1) + 1) = u;
    end
    
    [curves, curveIds, eidx] = traverse(vtx1, lst);
    
    totCurves = length(curves);
    if outputObj
        fout = fopen("slipstitchrib.obj", "wt");
        for i = 1 : size(vtx1, 1)
            fprintf(fout, "v %.6f %.6f %.6f\n", vtx1(i, 1), vtx1(i, 2), vtx1(i, 3));
        end
        for i = 1 : totCurves
            m = size(curves{i}, 1);
            fprintf(fout, "l");
            for j = 1 : m
                fprintf(fout, " %d", curveIds{i}{j});
            end
            fprintf(fout, "\n");
        end
        fclose(fout);
    end
      
    if outputMts
        fout = fopen("slipstitchrib.txt", "wt");
        fprintf(fout, "%d\n\n", totCurves);
        for i = 1 : totCurves
            m = size(curves{i}, 1);
            fprintf(fout, "%d\n", m);
            for j = 1 : m
                fprintf(fout, "%.6f %.6f %.6f\n", curves{i}(j, 1), curves{i}(j, 2), curves{i}(j, 3));
            end
            if i < totCurves
                fprintf(fout, "\n");
            end
        end
        fclose(fout);
    end
    
    if plotResult
        clf
        hold on
        for i = 1 : totCurves
            plot3(curves{i}(:, 1), curves{i}(:, 2), curves{i}(:, 3))
        end
        scatter3(vtx1(eidx, 1), vtx1(eidx, 2), vtx1(eidx, 3))
        hold off
        axis equal
    end
end


function [curves, curveIds, eidx] = traverse(vtx1, lst)
    totVtx1 = size(vtx1, 1);

    curves = {};
    curveIds = {};
    totCurves = 0;
    used = false(1, totVtx1);
    eidx = [];
    for i = 1 : totVtx1
        if lst(i, 1) == 1 && ~used(i)
            curCurve = {};
            j = i; t = 0;
            while true
                used(j) = true;
                t = t + 1;
                curCurve{t} = j; %#ok<AGROW>
                j1 = 0;
                for k = 1 : lst(j, 1)
                    if ~used(lst(j, k + 1))
                        j1 = lst(j, k + 1);
                        break;
                    end
                end
                
                if j1 > 0
                    j = j1;
                else
                    break;
                end
            end
            
            totCurves = totCurves + 1;
            curves{totCurves} = zeros(t, 3); %#ok<AGROW>
            for j = 1 : t
                curves{totCurves}(j, :) = vtx1(curCurve{j}, :);
            end
            curveIds{totCurves} = curCurve; %#ok<AGROW>
            eidx = [eidx, curCurve{1}, curCurve{t}]; %#ok<AGROW>
        end
    end
    assert(all(used))
end

