x = read_csv("E:\College\Courseware\Bachelor Thesis\programme\PCA2d plot\fvx.txt");
y = read_csv("E:\College\Courseware\Bachelor Thesis\programme\PCA2d plot\fvy.txt");
a = read_csv("E:\College\Courseware\Bachelor Thesis\programme\PCA2d plot\cbx.txt");
b = read_csv("E:\College\Courseware\Bachelor Thesis\programme\PCA2d plot\cby.txt");

x = evstr(x);
y = evstr(y);
a = evstr(a);
b = evstr(b);


plot([x],y,".r");

curr_entity = gce();
curr_polyline = curr_entity.children;
curr_polyline.mark_size = 2;

plot([a],b,".k");

legend("fv", "cb");
