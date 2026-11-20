import pdfplumber, numpy as np, json, collections, sys
pdf=pdfplumber.open(sys.argv[1]); pg=pdf.pages[7]
cx=lambda w:(w['x0']+w['x1'])/2; cy=lambda w:(w['top']+w['bottom'])/2
W=[w for w in pg.extract_words() if w['top']>200]
xl=[w for w in W if w['text'] in ('0','500','1000','1500') and (470<w['top']<490 or 678<w['top']<695)]
panels=[]
for row,(lo,hi) in enumerate([(470,490),(678,695)]):
    r=sorted([w for w in xl if lo<w['top']<hi],key=cx)
    for k in range(0,len(r),4):
        g=r[k:k+4]; v=[float(a['text']) for a in g]; assert v==[0,500,1000,1500],v
        panels.append(dict(row=row,xc=np.polyfit([cx(a) for a in g],v,1),x0=cx(g[0]),x1=cx(g[3])))
yl=[w for w in W if w['text'] in ('0','500','1000','1500','2000','2500','3000','3500')]
for p in panels:
    c=[w for w in yl if p['x0']-25<w['x1']<p['x0']-3 and ((p['row']==0 and 300<w['top']<475) or (p['row']==1 and 505<w['top']<680))]
    c=sorted(c,key=cy); p['yc']=np.polyfit([cy(a) for a in c],[float(a['text']) for a in c],1)
cs=[c for c in pg.curves if c['top']>250 and c.get('fill') and (c['x1']-c['x0'])<3.5]
ded=[]
for c in cs:
    mx=(c['x0']+c['x1'])/2; my=(c['top']+c['bottom'])/2
    sex='M' if tuple(c['non_stroking_color'])[0]<0.5 else 'F'
    row=0 if my<480 else 1
    pi=[i for i,q in enumerate(panels) if q['row']==row and q['x0']-5<=mx<=q['x1']+25][0]
    if any(d['panel']==pi and d['sex']==sex and abs(d['mx']-mx)<0.3 and abs(d['my']-my)<0.3 for d in ded): continue
    ded.append(dict(panel=pi,sex=sex,P=float(np.polyval(panels[pi]['xc'],mx)),C=float(np.polyval(panels[pi]['yc'],my)),mx=mx,my=my))
print(collections.Counter((o['panel'],o['sex']) for o in ded))
L=[l for l in pg.lines if l['top']>250]
res=[]
for pi in range(6):
    for s in 'MF':
        g=[o for o in ded if o['panel']==pi and o['sex']==s]; e=max(g,key=lambda o:o['C']); p=panels[pi]
        v=[l for l in L if abs(l['x0']-l['x1'])<0.2 and abs(l['x0']-e['mx'])<0.5 and l['top']-1<e['my']<l['bottom']+1 and (l['bottom']-l['top'])>0.5]
        h=[l for l in L if abs(l['top']-l['bottom'])<0.2 and abs(l['top']-e['my'])<0.5 and l['x0']-1<e['mx']<l['x1']+1 and (l['x1']-l['x0'])>0.5]
        cse=(max(l['bottom'] for l in v)-min(l['top'] for l in v))/2*abs(p['yc'][0]) if v else float('nan')
        pse=(max(l['x1'] for l in h)-min(l['x0'] for l in h))/2*abs(p['xc'][0]) if h else float('nan')
        res.append(dict(panel=pi+1,sex=s,n_points=len(g),P_ug=round(e['P'],1),C_ug=round(e['C'],1),P_se=round(pse,1),C_se=round(cse,1)))
for r in res: print(r)
json.dump(res,open('jensen_endpoints.json','w'))
