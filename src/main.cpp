#include "main.h"
#include "node.h"
#include "scene.h"
#include "SceneOcclusion.h"
#include "kernel.h"
#include "KernelSampler.h"
#include <python3.5m/Python.h>

#include <unordered_map>
#include <functional> //std::function
#include <thread> //std::thread

#include <time.h> //logging gimmick
#include <stdarg.h>

static RenderKernel *gpkernel = 0;
static Scene *gpscene = 0;
static SceneOcclusion *gpsceneocc = 0;
static enum ENGINE_STATE{
	ENGINE_STATE_READY,
	ENGINE_STATE_PROCESSING
} gstate = ENGINE_STATE_READY; //for asynchronous Python-scripting, to prevent locking up the Blender interface

#ifndef __unix__
#define strcasecmp stricmp
#endif

void DebugPrintf(const char *pfmt, ...){
	time_t rt;
	time(&rt);
	const struct tm *pti = localtime(&rt);

	char tbuf[256];
	strftime(tbuf,sizeof(tbuf),"[Droplet %F %T]",pti);
	printf("%s ",tbuf);

	va_list args;
	va_start(args,pfmt);
	vprintf(pfmt,args);
	va_end(args);
}

inline float PyGetFloat(PyObject *pb, const char *pn){
	PyObject *pa = PyObject_GetAttrString(pb,pn);
	float r = (float)PyFloat_AsDouble(pa);
	Py_DECREF(pa);
	return r;
}

inline uint PyGetUint(PyObject *pb, const char *pn){
	PyObject *pa = PyObject_GetAttrString(pb,pn);
	uint r = (uint)PyLong_AsLong(pa);
	Py_DECREF(pa);
	return r;
}

static PyObject * DRE_BeginRender(PyObject *pself, PyObject *pargs){
	PyObject *pscene, *pdata;
	uint tilex, tiley, w, h, smask;
	if(gpscene || !PyArg_ParseTuple(pargs,"OOIIIII",&pscene,&pdata,&tilex,&tiley,&w,&h,&smask)){
		DebugPrintf("Invalid arguments\n");
		return 0;
	}
	PyObject *pycam = PyObject_GetAttrString(pscene,"camera");

	PyObject *pyloc = PyObject_GetAttrString(pycam,"location");
	const float4 u = float4(0,0,-1,0);
	float4 p = float4(
		PyGetFloat(pyloc,"x"),
		PyGetFloat(pyloc,"y"),
		PyGetFloat(pyloc,"z"),1.0f);
	PyObject *pyrot = PyObject_GetAttrString(pycam,"rotation_euler"); //rotation_quaternion
	PyObject *pyquat = PyObject_CallMethod(pyrot,"to_quaternion","");
	float4 q = float4(
		PyGetFloat(pyquat,"x"),
		PyGetFloat(pyquat,"y"),
		PyGetFloat(pyquat,"z"),
		PyGetFloat(pyquat,"w"));
	float4 d = XMVector3Rotate(u.v,q.v);
	Py_DECREF(pyloc);
	Py_DECREF(pyrot);
	Py_DECREF(pyquat);

	PyObject *pycdt = PyObject_GetAttrString(pycam,"data");

	PyObject *pyvfrm = PyObject_CallMethod(pycdt,"view_frame","O",pscene);
	PyObject *pyvfrmv = PyObject_GetIter(pyvfrm);
	PyObject *pyvfrma = PyIter_Next(pyvfrmv);
	float frx = PyGetFloat(pyvfrma,"x");
	float fry = PyGetFloat(pyvfrma,"y");
	Py_DECREF(pyvfrma);
	Py_DECREF(pyvfrmv);
	Py_DECREF(pyvfrm);

	//fov = 2*atanf((0.5*sensor_size)/bcam->lens/aspectratio) = camera.data.angle;
	//float fov = 1.045399f*(float)h/(float)w*PyGetFloat(pycdt,"angle");
	float zmin = PyGetFloat(pycdt,"clip_start");
	float zmax = PyGetFloat(pycdt,"clip_end");

	Py_DECREF(pycdt);
	Py_DECREF(pycam);

	XMMATRIX view = XMMatrixLookToRH(p.v,d.v,u.v);
	XMMATRIX proj = XMMatrixPerspectiveRH(frx,fry,zmin,zmax); //XMMatrixPerspectiveFovRH(fov,(float)w/(float)h,zmin,zmax);

	dmatrix44 sviewi, sproji;
	matrix44::store(&sviewi,matrix44(XMMatrixInverse(0,view).r));
	matrix44::store(&sproji,matrix44(XMMatrixInverse(0,proj).r));

	std::unordered_map<Py_hash_t, Node::NodeTree *> ntm;

	PyObject *pngn = PyObject_GetAttrString(pdata,"node_groups");
	PyObject *pntl = PyObject_CallMethod(pngn,"values",""); //TODO: iterator
	uint ntc = PyList_Size(pntl);
	for(uint i = 0; i < ntc; ++i){
		PyObject *pnt1 = PyList_GetItem(pntl,i);
		PyObject *pns1 = PyObject_GetAttrString(pnt1,"nodes");

		std::unordered_map<Py_hash_t, Node::BaseNode *> nodem; //keep list of local duplicates
		std::function<Node::BaseNode * (PyObject *, uint, Node::NodeTree *)> ntree = [&](PyObject *proot, uint l, Node::NodeTree *pnt)->Node::BaseNode *{
			//
			PyObject *pidn = PyObject_GetAttrString(proot,"bl_idname"); //TODO: give enum id
			const char *pn = PyUnicode_AsUTF8(pidn);
			//printf("%s\n",pn);

			Node::BaseNode *pbn = Node::CreateNodeByType(pn,proot,l,pnt);//Node::CreateNode(pn,l);
			if(!pbn){
				DebugPrintf("Error: unknown node %s\n",pn);
				return 0;
			}

			Py_DECREF(pidn);

			PyObject *pnouts = PyObject_GetAttrString(proot,"outputs");
			PyObject *pnoutv = PyObject_GetIter(pnouts);
			//
			uint nx = 0;
			for(PyObject *pnout1 = PyIter_Next(pnoutv); pnout1; Py_DecRef(pnout1), pnout1 = PyIter_Next(pnoutv), ++nx){
				//
				PyObject *plnks = PyObject_GetAttrString(pnout1,"links");
				if(PyTuple_Size(plnks) > 0)
					pbn->omask |= 1<<nx;
				Py_DECREF(plnks);
			}

			Py_DECREF(pnoutv);
			Py_DECREF(pnouts);

			Py_hash_t h = PyObject_Hash(proot);
			nodem.insert(std::pair<Py_hash_t, Node::BaseNode *>(h,pbn));

			PyObject *pnins = PyObject_GetAttrString(proot,"inputs");
			PyObject *pninv = PyObject_GetIter(pnins);
			//
			nx = 0;
			for(PyObject *pnin1 = PyIter_Next(pninv); pnin1; Py_DecRef(pnin1), pnin1 = PyIter_Next(pninv), ++nx){
				PyObject *plnks = PyObject_GetAttrString(pnin1,"links");
				PyObject *plnkv = PyObject_GetIter(plnks);
				PyObject *plnk1 = PyIter_Next(plnkv);
				Py_DECREF(plnks);
				Py_DECREF(plnkv);

				pidn = PyObject_GetAttrString(pnin1,"bl_idname");
				pn = PyUnicode_AsUTF8(pidn);

				if(!plnk1){
					PyObject *pvalue = PyObject_GetAttrString(pnin1,"value"); //note: every socket class should have this property, even if not used
					pbn->pnodes[nx] = Node::CreateNodeBySocket(pn,pvalue,l+1,pnt);
					pbn->indices[nx] = 0;

					Py_DECREF(pvalue);
					Py_DECREF(pidn);
					continue;
				}else pbn->imask |= 1<<nx;

				PyObject *pnode = PyObject_GetAttrString(plnk1,"from_node");
				Py_hash_t h = PyObject_Hash(pnode);

				PyObject *psock = PyObject_GetAttrString(plnk1,"from_socket");
				Py_hash_t t = PyObject_Hash(psock);

				//get the child node output index
				PyObject *pcouts = PyObject_GetAttrString(pnode,"outputs");
				PyObject *pcoutv = PyObject_GetIter(pcouts);
				uint sx = 0;
				for(PyObject *pcout1 = PyIter_Next(pcoutv); pcout1; Py_DecRef(pcout1), pcout1 = PyIter_Next(pcoutv)){
					Py_hash_t q = PyObject_Hash(pcout1);
					if(t == q)
						break;
					PyObject *pidn1 = PyObject_GetAttrString(pcout1,"bl_idname");
					const char *pn1 = PyUnicode_AsUTF8(pidn1);
					if(strcmp(pn,pn1) == 0)
						++sx;
				}
				pbn->indices[nx] = sx;

				Py_DECREF(pcoutv);
				Py_DECREF(pcouts);
				//
				Py_DECREF(psock);

				std::unordered_map<Py_hash_t, Node::BaseNode *>::const_iterator m = nodem.find(h);
				if(m != nodem.end()){
					if(m->second->level < l+1)
						m->second->level = l+1;
					pbn->pnodes[nx] = m->second;
				}else pbn->pnodes[nx] = ntree(pnode,l+1,pnt);

				Py_DECREF(plnk1);
				Py_DECREF(pnode);
			}

			Py_DECREF(pninv);
			Py_DECREF(pnins);

			return pbn;
		};

		PyObject *pidn = PyObject_GetAttrString(pnt1,"bl_idname"); //TODO: give enum id
		const char *pn = PyUnicode_AsUTF8(pidn);
		if(strcmp(pn,"ClNodeTree") != 0){
			Py_DECREF(pidn);
			Py_DECREF(pns1);
			continue;
		}
		Py_DECREF(pidn);

		PyObject *pno = PyObject_GetAttrString(pnt1,"name");
		Node::NodeTree *pntree = new Node::NodeTree(PyUnicode_AsUTF8(pno));
		Py_DECREF(pno);

		PyObject *proot = PyObject_CallMethod(pns1,"get","(s)","Surface Output");
		if(!proot){
			DebugPrintf("Warning: output node not found\n");
			//Py_DECREF(proot);
			Py_DECREF(pns1);
			continue;
		}
		ntree(proot,0,pntree);

		pntree->SortNodes();
		pntree->ApplyBranchMask();

		Py_hash_t h = PyObject_Hash(pnt1);
		ntm.insert(std::pair<Py_hash_t, Node::NodeTree *>(h,pntree));

		Py_DECREF(pns1);
		//Py_DECREF(pnt1); //borrowed ref
	}

	Py_DECREF(pngn);

	PyObject *pbmeshn = PyUnicode_FromString("bmesh");
	PyObject *pbmesh = PyImport_Import(pbmeshn);

	PyObject *pbmops = PyObject_GetAttrString(pbmesh,"ops");

	PyObject *pyoat = PyObject_GetAttrString(pscene,"objects");
	PyObject *pyobl = PyObject_CallMethod(pyoat,"values","");
	uint objc = PyList_Size(pyobl);
	for(uint i = 0; i < objc; ++i){
		PyObject *pobj = PyList_GetItem(pyobl,i);

		//check whether the object is on an active render layer
		PyObject *pyrls = PyObject_GetAttrString(pobj,"layers");
		PyObject *pyrlm = PyObject_GetIter(pyrls);
		uint rx = 0, omask = 0;
		for(PyObject *pyrl1 = PyIter_Next(pyrlm); pyrl1; Py_DecRef(pyrl1), pyrl1 = PyIter_Next(pyrlm), ++rx)
			omask |= PyObject_IsTrue(pyrl1)<<rx;
		Py_DECREF(pyrlm);
		Py_DECREF(pyrls);
		if((smask&omask) == 0)
			continue;

		//get particle systems
		PyObject *ppro = PyObject_GetAttrString(pobj,"particle_systems");
		PyObject *ppsl = PyObject_CallMethod(ppro,"values","");
		uint prc = PyList_Size(ppsl);
		for(uint j = 0; j < prc; ++j){
			PyObject *pps = PyList_GetItem(ppsl,j);
			PyObject *pst = PyObject_GetAttrString(pps,"settings"); //settings.droplet?
			PyObject *psd = PyObject_GetAttrString(pst,"droplet");
			PyObject *ppn = PyObject_GetAttrString(psd,"nodetree");
			const char *pname = PyUnicode_AsUTF8(ppn);

			std::unordered_map<Py_hash_t, Node::NodeTree *>::const_iterator m = std::find_if(ntm.begin(),ntm.end(),[=](const std::unordered_map<Py_hash_t, Node::NodeTree *>::value_type &t)->bool{
				return strcmp(t.second->name,pname) == 0;
			});
			if(m == ntm.end()){
				DebugPrintf("Warning: invalid node tree %s. Skipping particle system (index=%u,%u).\n",pname,i,j);
				continue;
			}
			SceneData::ParticleSystem *pprs = new SceneData::ParticleSystem(m->second); //TODO: reserve() the particle vector size

			Py_DECREF(ppn);
			Py_DECREF(psd);
			Py_DECREF(pst);

			//TODO: Get visibility status. There's some weird modifier class for this.
			PyObject *ppr = PyObject_GetAttrString(pps,"particles");
			PyObject *ppi = PyObject_GetIter(ppr);
			for(PyObject *pni = PyIter_Next(ppi); pni; Py_DecRef(pni), pni = PyIter_Next(ppi)){
				PyObject *past = PyObject_GetAttrString(pni,"alive_state");
				const char *past1 = PyUnicode_AsUTF8(past);
				if(strcasecmp(past1,"DEAD") == 0){
					Py_DECREF(past);
					continue;
				}
				Py_DECREF(past);

				PyObject *pvco = PyObject_GetAttrString(pni,"location");
				float4 co = float4(
					PyGetFloat(pvco,"x"),
					PyGetFloat(pvco,"y"),
					PyGetFloat(pvco,"z"),1.0f);
				pprs->pl.push_back(dfloat3(co)); //already in world-space
				Py_DECREF(pvco);

				PyObject *pvvl = PyObject_GetAttrString(pni,"velocity");
				float4 vl = float4(
					PyGetFloat(pvvl,"x"),
					PyGetFloat(pvvl,"y"),
					PyGetFloat(pvvl,"z"),1.0f);
				pprs->vl.push_back(dfloat3(vl));
				Py_DECREF(pvvl);
			}
			Py_DECREF(ppi);
			Py_DECREF(ppr);
		}
		Py_DECREF(ppro);

		if(PyGetUint(pobj,"hide_render") != 0)
			continue;

		PyObject *ptype = PyObject_GetAttrString(pobj,"type");
		const char *ptype1 = PyUnicode_AsUTF8(ptype);

		if(strcasecmp(ptype1,"LAMP") == 0){
			pyrot = PyObject_GetAttrString(pobj,"rotation_euler");
			pyquat = PyObject_CallMethod(pyrot,"to_quaternion","");

			q = float4(
				PyGetFloat(pyquat,"x"),
				PyGetFloat(pyquat,"y"),
				PyGetFloat(pyquat,"z"),
				PyGetFloat(pyquat,"w"));
			d = float4(XMVector3Rotate(u.v,q.v));

			Py_DECREF(pyrot);
			Py_DECREF(pyquat);

			PyObject *pdata = PyObject_GetAttrString(pobj,"data");
			PyObject *psd = PyObject_GetAttrString(pdata,"droplet");
			PyObject *pcolor = PyObject_GetAttrString(psd,"color");
			float intensity = PyGetFloat(psd,"intensity");

			dfloat3 D = dfloat3(-d);
			dfloat3 c = intensity*dfloat3(
				PyGetFloat(pcolor,"r"),
				PyGetFloat(pcolor,"g"),
				PyGetFloat(pcolor,"b"));
			float angle = PyGetFloat(psd,"angle");
			KernelSampler::BaseLight *pl = new KernelSampler::SunLight(&D,&c,angle);

			Py_DECREF(pcolor);
			Py_DECREF(psd);
			Py_DECREF(pdata);
		}

		else if(strcasecmp(ptype1,"MESH") == 0){
			PyObject *psd = PyObject_GetAttrString(pobj,"droplet");
			PyObject *ppn = PyObject_GetAttrString(psd,"nodetree");
			PyObject *pyholdout = PyObject_GetAttrString(psd,"holdout");
			const char *pname = PyUnicode_AsUTF8(ppn);
			bool holdout = PyObject_IsTrue(pyholdout);

			std::unordered_map<Py_hash_t, Node::NodeTree *>::const_iterator m = std::find_if(ntm.begin(),ntm.end(),[=](const std::unordered_map<Py_hash_t, Node::NodeTree *>::value_type &t)->bool{
				return strcmp(t.second->name,pname) == 0;
			});
			if(m == ntm.end()){
				DebugPrintf("Warning: invalid node tree %s. Skipping surface (index=%u).\n",pname,i);
				continue;
			}
			SceneData::Surface *psobj = new SceneData::Surface(m->second,holdout);

			//check if SmokeCache node was used
			for(uint j = 0; j < m->second->nodes1.size(); ++j){
				Node::ISmokeCache *pscn = dynamic_cast<Node::ISmokeCache *>(m->second->nodes1[j]);
				if(pscn){
					PyObject *pyvdb = PyObject_GetAttrString(psd,"vdbcache");
					PyObject *pyrho = PyObject_GetAttrString(psd,"vdbrho");
					PyObject *pyvel = PyObject_GetAttrString(psd,"vdbvel");
					SceneData::SmokeCache *pprs = new SceneData::SmokeCache(m->second,
						PyUnicode_AsUTF8(pyvdb),PyUnicode_AsUTF8(pyrho), PyUnicode_AsUTF8(pyvel));
					Py_DECREF(pyvdb);
					Py_DECREF(pyrho);
					Py_DECREF(pyvel);
					break;
				}
			}

			Py_DECREF(pyholdout);
			Py_DECREF(ppn);
			Py_DECREF(psd);

			XMMATRIX wm;
			PyObject *pwm = PyObject_GetAttrString(pobj,"matrix_world");
			PyObject *pmr = PyObject_GetAttrString(pwm,"row");
			PyObject *pmi = PyObject_GetIter(pmr);
			uint mr = 0;
			for(PyObject *pni = PyIter_Next(pmi); pni; Py_DecRef(pni), pni = PyIter_Next(pmi), ++mr){
				wm.r[mr] = float4(
					PyGetFloat(pni,"x"),
					PyGetFloat(pni,"y"),
					PyGetFloat(pni,"z"),
					PyGetFloat(pni,"w")).v;
			}
			Py_DECREF(pmi);
			Py_DECREF(pmr);
			Py_DECREF(pwm);
			wm = XMMatrixTranspose(wm);

			PyObject *pbm = PyObject_CallMethod(pbmesh,"new","");
			PyObject_CallMethod(pbm,"from_object","OO",pobj,pscene);

			PyObject *pyfa = PyObject_GetAttrString(pbm,"faces");
			PyObject *pyag = Py_BuildValue("(O)",pbm);
			PyObject *pykw = PyDict_New();
			PyObject *pytm = PyObject_GetAttrString(pbmops,"triangulate");
			PyDict_SetItemString(pykw,"faces",pyfa);
			PyObject_Call(pytm,pyag,pykw);
			Py_DECREF(pyfa);
			Py_DECREF(pyag);
			Py_DECREF(pytm);
			Py_DECREF(pykw);

			PyObject *pvs = PyObject_GetAttrString(pbm,"verts");
			PyObject *pvi = PyObject_GetIter(pvs);
			for(PyObject *pni = PyIter_Next(pvi); pni; Py_DecRef(pni), pni = PyIter_Next(pvi)){
				PyObject *pvco = PyObject_GetAttrString(pni,"co");
				float4 co = float4(
					PyGetFloat(pvco,"x"),
					PyGetFloat(pvco,"y"),
					PyGetFloat(pvco,"z"),1.0f);
				dfloat3 sco = dfloat3(XMVector3TransformCoord(co.v,wm));
				psobj->vl.push_back(sco);

				Py_DECREF(pvco);
			}

			Py_DECREF(pvi);
			Py_DECREF(pvs);

			PyObject *pfs = PyObject_GetAttrString(pbm,"faces");
			PyObject *pfi = PyObject_GetIter(pfs);
			for(PyObject *pni = PyIter_Next(pfi); pni; Py_DecRef(pni), pni = PyIter_Next(pfi)){
				PyObject *pfvs = PyObject_GetAttrString(pni,"verts");
				PyObject *pfvi = PyObject_GetIter(pfvs);
				for(PyObject *pji = PyIter_Next(pfvi); pji; Py_DecRef(pji), pji = PyIter_Next(pfvi)){
					uint index = PyGetUint(pji,"index");
					psobj->tl.push_back(index);
				}
				Py_DECREF(pfvi);
			}

			Py_DECREF(pfi);
			Py_DECREF(pfs);

			PyObject_CallMethod(pbm,"free","");
			Py_DECREF(pbm);

			//Py_DECREF(pobj); //borrowed ref
		}

		Py_DECREF(ptype);
	}

	Py_DECREF(pyoat);

	Py_DECREF(pbmops);
	Py_DECREF(pbmesh);

	PyObject *pysampling = PyObject_GetAttrString(pscene,"blcloudsampling");
	uint scattevs = PyGetUint(pysampling,"scatterevs");
	float msigmas = PyGetFloat(pysampling,"msigmas");
	float msigmaa = PyGetFloat(pysampling,"msigmaa");
	PyObject *pypf = PyObject_GetAttrString(pysampling,"phasef");
	const char *pypfs = PyUnicode_AsUTF8(pypf);
	KernelSampler::PhaseFunction *ppf;
	switch(pypfs[0]){
	case 'H':
		ppf = &KernelSampler::HGPhase::ghg;
		DebugPrintf("Using Henyey-Greenstein phase.\n");
		break;
	case 'M':
		ppf = &KernelSampler::MiePhase::gmie;
		DebugPrintf("Using Mie phase.\n");
		break;
	}
	Py_DECREF(pypf);
	Py_DECREF(pysampling);

	PyObject *pygrid = PyObject_GetAttrString(pscene,"blcloudgrid");
	float dsize = PyGetFloat(pygrid,"detailsize");
	float qband = PyGetFloat(pygrid,"qfbandw");
	uint maxd = PyGetUint(pygrid,"maxdepth");
	Py_DECREF(pygrid);

	PyObject *pyperf = PyObject_GetAttrString(pscene,"blcloudperf");
	PyObject *pycache = PyObject_GetAttrString(pyperf,"cache");
	const char *pcms = PyUnicode_AsUTF8(pycache);
	SCENE_CACHE_MODE cm = (pcms[0] == 'R'?SCENE_CACHE_READ:pcms[0] == 'W'?SCENE_CACHE_WRITE:SCENE_CACHE_DISABLED);

	Py_DECREF(pycache);
	Py_DECREF(pyperf);

	//TODO: cache postfix from bpy.path.basename(bpy.data.filepath)

	PyObject *pyworld = PyObject_GetAttrString(pscene,"world");
	PyObject *pywsettings = PyObject_GetAttrString(pyworld,"droplet");
	PyObject *pyocclusion = PyObject_GetAttrString(pywsettings,"occlusion");
	bool occlusion = PyObject_IsTrue(pyocclusion);

	Py_DECREF(pyocclusion);
	Py_DECREF(pywsettings);
	Py_DECREF(pyworld);

	gstate = ENGINE_STATE_PROCESSING;

	std::thread async([=]()->void{
		if(occlusion){
			gpsceneocc = new SceneOcclusion();
			gpsceneocc->Initialize();
		}

		gpscene = new Scene(); //TODO: interface for blender status reporting
		gpscene->Initialize(dsize,maxd,qband,cm);

		gpkernel = new RenderKernel();
		gpkernel->Initialize(gpscene,gpsceneocc,&sviewi,&sproji,ppf,scattevs,msigmas,msigmaa,tilex,tiley,w,h,0);

		SceneData::SmokeCache::DeleteAll();
		SceneData::ParticleSystem::DeleteAll();
		SceneData::Surface::DeleteAll();
		Node::NodeTree::DeleteAll();

		gstate = ENGINE_STATE_READY;
	});
	async.detach();

	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject * DRE_Render(PyObject *pself, PyObject *pargs){
	uint x0, y0, tilex, tiley, samples;
	if(!gpkernel || !PyArg_ParseTuple(pargs,"IIIII",&x0,&y0,&tilex,&tiley,&samples))
		return 0;

	Py_INCREF(Py_None);
	if(gstate != ENGINE_STATE_READY)
		return Py_None;

	gstate = ENGINE_STATE_PROCESSING;
	std::thread async([=]()->void{
		gpkernel->Render(x0,y0,tilex,tiley,samples);
		gstate = ENGINE_STATE_READY;
	});
	async.detach();

	return Py_None;
}

static PyObject * DRE_EndRender(PyObject *pself, PyObject *pargs){
	KernelSampler::BaseLight::DeleteAll();

	gpkernel->Destroy();
	delete gpkernel;
	gpkernel = 0;

	gpscene->Destroy();
	delete gpscene;
	gpscene = 0;

	if(gpsceneocc){
		gpsceneocc->Destroy();
		delete gpsceneocc;
		gpsceneocc = 0;
	}

	DebugPrintf("Shutdown\n");

	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject * DRE_QueryStatus(PyObject *pself, PyObject *pargs){
	//
	return Py_BuildValue("i",gstate);
}

static PyObject * DRE_QueryResult(PyObject *pself, PyObject *pargs){
	if(!gpkernel)
		return 0;
	if(gstate != ENGINE_STATE_READY){
		Py_INCREF(Py_None);
		return Py_None;
	}

	uint l = gpkernel->tilew*gpkernel->tileh;
	PyObject *prt = PyList_New(l);
	for(uint i = 0; i < l; ++i){
		PyObject *pc = Py_BuildValue("[f,f,f,f]",gpkernel->phb[i].x,gpkernel->phb[i].y,gpkernel->phb[i].z,gpkernel->phb[i].w);
		PyList_SET_ITEM(prt,i,pc);
	}

	return prt;
}

static PyMethodDef g_blmethods[] = {
	{"BeginRender",DRE_BeginRender,METH_VARARGS,"BeginRender() doc string"}, //CreateDevice
	{"Render",DRE_Render,METH_VARARGS,"Render() doc string"},
	{"EndRender",DRE_EndRender,METH_NOARGS,"EndRender() doc string"},
	{"QueryStatus",DRE_QueryStatus,METH_NOARGS,"QueryStatus() doc string"},
	{"QueryResult",DRE_QueryResult,METH_NOARGS,"QueryResult() doc string"},
	{0,0,0,0}
};

static struct PyModuleDef g_bldre = {
	PyModuleDef_HEAD_INIT,
	"droplet",
	"Droplet Render Engine",
	-1,
	g_blmethods,
	0,0,0,0 //freefunc m_free
};

#define BLCLOUD_MODINITFUNC PyInit_libdroplet

PyMODINIT_FUNC BLCLOUD_MODINITFUNC(){
	return PyModule_Create(&g_bldre);
}
