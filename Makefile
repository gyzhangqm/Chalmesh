include makefile.defs

chalmesh:
	cd $(UTIL);     $(MAKE); cd ..
	cd $(OGL_PLOT); $(MAKE); cd ..
	cd $(STRETCH);  $(MAKE); cd ..
	cd $(SURFACE);  $(MAKE); cd ..
	cd $(HOLE);     $(MAKE); cd ..
	cd $(VOLUME);   $(MAKE); cd ..
	cd $(OVERLAP);  $(MAKE); cd ..
	if [ "$(MESA_HOME)" = "" ]; then \
          cd $(MAIN); $(MAKE) ../bin/chalmesh.$(SYSTEM); cd ..; \
        else \
          cd $(MAIN); $(MAKE) ../bin/chalmesh_mesa.$(SYSTEM); cd ..; \
        fi
	cd $(DATA);     $(MAKE); cd ..

clean:
	cd $(UTIL);     $(MAKE) clean; cd ..
	cd $(OGL_PLOT); $(MAKE) clean; cd ..
	cd $(STRETCH);  $(MAKE) clean; cd ..
	cd $(SURFACE);  $(MAKE) clean; cd ..
	cd $(HOLE);     $(MAKE) clean; cd ..
	cd $(VOLUME);   $(MAKE) clean; cd ..
	cd $(OVERLAP);  $(MAKE) clean; cd ..
	cd $(MAIN);     $(MAKE) clean; cd ..


