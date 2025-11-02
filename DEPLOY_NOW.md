# 🚀 Deploy CriOS to Railway RIGHT NOW

## You're Ready! All Configuration Complete ✅

Everything is configured and committed to your branch: `claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`

---

## 🎯 Deploy in 3 Ways

### Option 1: Automated Script (Fastest - 5 minutes)

Run this from your **local machine** where you have the repository:

```bash
# 1. Pull the latest changes
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
git pull origin claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH

# 2. Install Railway CLI
npm i -g @railway/cli

# 3. Run the deployment script
./deploy-railway.sh
```

The script will:
- ✅ Login to Railway
- ✅ Create project
- ✅ Add PostgreSQL & Redis
- ✅ Deploy all 6 services
- ✅ Run health checks

---

### Option 2: Railway Dashboard (Easiest - 2 clicks)

**Perfect if you prefer UI over CLI:**

1. **Go to Railway**: https://railway.app
2. **Click "New Project"**
3. **Select "Deploy from GitHub repo"**
4. **Choose**: `MichaelCrowe11/crios-dr-crowe-coder`
5. **Select branch**: `claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`
6. **Railway auto-detects and deploys!** 🎉

Then add databases:
- Click "+ New" → "Database" → "PostgreSQL"
- Click "+ New" → "Database" → "Redis"

Set required environment variables in each service:
- `ANTHROPIC_API_KEY=your_key_here`
- `JWT_SECRET` (Railway can generate this)

**That's it!** Railway handles everything else automatically.

---

### Option 3: Manual CLI (Most Control - 10 minutes)

From your **local machine**:

```bash
# 1. Pull latest changes
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
git pull

# 2. Install Railway CLI
npm i -g @railway/cli

# 3. Login
railway login

# 4. Create project
railway init

# Give it a name like: crios-dr-crowe-coder

# 5. Add PostgreSQL
railway add --plugin postgresql

# 6. Add Redis
railway add --plugin redis

# 7. Set critical environment variables
railway variables set ANTHROPIC_API_KEY=sk-ant-your_key_here
railway variables set JWT_SECRET=$(openssl rand -base64 32)
railway variables set ENVIRONMENT=production
railway variables set DEBUG=false

# 8. Deploy all services
# Railway will auto-detect the configuration from railway.toml
railway up

# 9. Check status
railway status

# 10. View your app!
railway open
```

---

## 📋 Before You Deploy - Quick Checklist

- [ ] **Railway Account**: Sign up at https://railway.app (free)
- [ ] **Anthropic API Key**: Get from https://console.anthropic.com
- [ ] **Git Branch**: Pull `claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`
- [ ] **Railway CLI** (if using CLI): `npm i -g @railway/cli`

---

## 🔑 Required Environment Variables

You'll need to set these in Railway:

### Critical (Required)
```bash
ANTHROPIC_API_KEY=sk-ant-your_key_here  # Get from console.anthropic.com
JWT_SECRET=your_random_32_char_secret   # Generate: openssl rand -base64 32
```

### Auto-Set by Railway (Don't set these manually)
```bash
DATABASE_URL=${{Postgres.DATABASE_URL}}
REDIS_URL=${{Redis.REDIS_URL}}
```

### Optional (For Advanced Features)
```bash
# Stripe (for payments)
STRIPE_SECRET_KEY=sk_live_...
STRIPE_PUBLIC_KEY=pk_live_...

# SendGrid (for emails)
SENDGRID_API_KEY=SG....

# OAuth (Google, GitHub, Microsoft)
GOOGLE_CLIENT_ID=...
GOOGLE_CLIENT_SECRET=...
```

Full list in: `.env.railway.template`

---

## 🎬 What Happens When You Deploy?

Railway will automatically:

1. **Detect Services** (from `railway.toml`):
   - Backend API (FastAPI + RDKit)
   - Frontend (Next.js 14)
   - Celery Worker (ML Pipeline)
   - Core API (Discovery Engine)

2. **Build Docker Images**:
   - Multi-stage builds
   - Install RDKit, PyTorch, ML dependencies
   - Optimize for production

3. **Provision Infrastructure**:
   - PostgreSQL database
   - Redis cache
   - Internal networking

4. **Deploy Services**:
   - Health checks enabled
   - Auto-scaling configured
   - SSL certificates (automatic)

5. **Give You URLs**:
   - Frontend: `https://crios-frontend-xxx.railway.app`
   - Backend: `https://crios-backend-xxx.railway.app`
   - Custom domains supported!

---

## 📊 Deployment Timeline

| Step | Time | What's Happening |
|------|------|------------------|
| Project Init | 30s | Creating Railway project |
| Database Setup | 1-2 min | Provisioning PostgreSQL & Redis |
| Backend Build | 3-5 min | Installing RDKit, Python deps |
| Frontend Build | 2-3 min | Building Next.js app |
| Worker Build | 2-3 min | Installing ML libraries |
| First Deploy | **~10 min** | All services starting |

**Subsequent deploys**: ~3-5 minutes (cached layers)

---

## 🔍 After Deployment - Verify

### Check Service Status
```bash
railway status
```

You should see:
```
✓ crios-backend     Deployed
✓ crios-frontend    Deployed
✓ crios-worker      Deployed
✓ postgres          Running
✓ redis             Running
```

### Test Backend API
```bash
# Get your backend URL from railway status
curl https://your-backend-url.railway.app/health

# Should return: {"status":"healthy"}
```

### Test Frontend
```bash
# Open in browser
railway open crios-frontend

# Or visit the URL from railway status
```

### View Logs
```bash
# Backend logs
railway logs --service crios-backend --follow

# Frontend logs
railway logs --service crios-frontend --follow

# Worker logs
railway logs --service crios-worker --follow
```

---

## 🐛 Common Issues & Quick Fixes

### Issue: "Build failed - RDKit not found"
**Solution**: This shouldn't happen with our Dockerfile, but if it does:
```dockerfile
# Ensure Dockerfile has:
RUN apt-get install -y librdkit-dev libboost-all-dev
```

### Issue: "Frontend can't connect to backend"
**Solution**: Set `NEXT_PUBLIC_API_URL`:
```bash
railway variables set NEXT_PUBLIC_API_URL=https://your-backend-url.railway.app --service crios-frontend
railway restart --service crios-frontend
```

### Issue: "Database connection refused"
**Solution**: Verify DATABASE_URL is set:
```bash
railway variables get DATABASE_URL
# Should start with: postgresql://...
```

### Issue: "Celery worker not processing jobs"
**Solution**: Check Redis connection:
```bash
railway logs --service crios-worker --follow
railway variables get REDIS_URL --service crios-worker
```

---

## 💰 Cost Breakdown

### Hobby Plan (~$10-15/month)
- ✅ 6 services (Backend, Frontend, Worker, Core API, PostgreSQL, Redis)
- ✅ 500 hours/month execution time
- ✅ 100GB bandwidth
- ✅ Shared CPU/RAM
- ✅ Community support

### Pro Plan (~$25-35/month)
- ✅ Everything in Hobby
- ✅ Better performance (dedicated resources)
- ✅ Custom domains
- ✅ Team collaboration
- ✅ Priority support
- ✅ Automatic backups

**First $5 is free!** (Hobby plan includes $5 credit)

---

## 🎓 Next Steps After Deployment

### Immediate
1. ✅ **Test all endpoints** - Use Postman or curl
2. ✅ **Add custom domain** - `railway domain add yourdomain.com`
3. ✅ **Set up monitoring** - Add Sentry DSN

### Short Term
1. ✅ **Run database migrations** - `railway run python manage.py migrate`
2. ✅ **Load sample data** - Test with demo compounds
3. ✅ **Configure OAuth** - Google, GitHub login

### Long Term
1. ✅ **Set up CI/CD** - Auto-deploy on git push
2. ✅ **Enable auto-backups** - Railway Pro feature
3. ✅ **Scale services** - Add replicas as needed
4. ✅ **Add staging env** - Clone project for testing

---

## 📚 Documentation Reference

| Doc | Use Case |
|-----|----------|
| `QUICK_START_RAILWAY.md` | 5-minute deploy guide |
| `RAILWAY_DEPLOYMENT.md` | Complete 30-page guide |
| `DEPLOYMENT_SUMMARY.md` | Integration overview |
| `.env.railway.template` | All env variables |
| `deploy-railway.sh` | Automated script |

---

## 🆘 Need Help?

### Railway Support
- **Docs**: https://docs.railway.app
- **Discord**: https://discord.gg/railway (very active!)
- **Status**: https://status.railway.app

### CriOS Support
- **GitHub Issues**: https://github.com/MichaelCrowe11/crios-dr-crowe-coder/issues
- **Documentation**: See `/docs` folder

---

## 🎉 You're Ready!

**Everything is configured and ready to deploy. Choose your method above and let's go! 🚀**

### Recommended: Use Option 2 (Railway Dashboard)
It's the easiest for first-time deployment:
1. Go to https://railway.app
2. "New Project" → "Deploy from GitHub"
3. Select your repo and branch
4. Done! ✅

**Happy Deploying! Your drug discovery platform will be live in ~10 minutes! 💊🔬**
