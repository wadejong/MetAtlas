mongod --dbpath /scratch/users/bkrull/db --auth > /dev/null &
export FWAPP_AUTH_USERNAME=metatlas
export FWAPP_AUTH_PASSWORD=metabolite
lpad -l ~/.fireworks/metatlas.yaml webgui --host=0.0.0.0 --port=8000 > /dev/null &
lpad -l local_db.yaml webgui --host=0.0.0.0 --port=8001 > /dev/null &
