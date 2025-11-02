# 🚀 START HERE - Deploy CriOS to Railway

## ⚡ Deploy in 3 Commands (From Your Computer)

```bash
# 1. Get the code
git clone https://github.com/MichaelCrowe11/crios-dr-crowe-coder.git
cd crios-dr-crowe-coder
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH

# 2. Run the automated deployment
./deploy-railway-cli.sh

# 3. Follow the prompts - Done in ~15 minutes! ✅
```

---

## 📋 Before You Start (2 minutes)

You need:

1. **Railway Account** (free): https://railway.app
2. **Anthropic API Key**: https://console.anthropic.com
3. **Node.js 18+**: Run `node -v` to check

---

## 🎯 Deployment Steps

### Step 1: Open Terminal on Your Computer

**Mac/Linux:**
- Press `Cmd+Space` and type "Terminal"
- Or find it in Applications > Utilities

**Windows:**
- Press `Win+R`, type `cmd`, press Enter
- Or use PowerShell or Git Bash

### Step 2: Run These Commands

Copy and paste one at a time:

```bash
# Navigate to where you want the code
cd ~/Documents  # or wherever you like

# Clone the repository
git clone https://github.com/MichaelCrowe11/crios-dr-crowe-coder.git

# Enter the directory
cd crios-dr-crowe-coder

# Checkout the deployment branch
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH

# Make the script executable (Mac/Linux)
chmod +x deploy-railway-cli.sh

# Run the deployment script
./deploy-railway-cli.sh
```

**On Windows:** If `./deploy-railway-cli.sh` doesn't work, use:
```bash
bash deploy-railway-cli.sh
```

### Step 3: Follow the Prompts

The script will ask you:

1. **"Ready to login to Railway?"** → Press `y`
   - Browser opens → Login to Railway

2. **"Create new Railway project?"** → Press `y`
   - Enter a name like "crios-dr-crowe-coder"

3. **"Add PostgreSQL database?"** → Press `y`
   - Wait ~1 minute while it provisions

4. **"Add Redis cache?"** → Press `y`
   - Wait ~1 minute while it provisions

5. **"Enter your Anthropic API key"** → Paste your key
   - Get it from: https://console.anthropic.com

6. **"Start deployment now?"** → Press `y`
   - Wait ~10-15 minutes (grab coffee ☕)

### Step 4: Done!

The script will show you:
- ✅ Your service URLs
- ✅ Health check results
- ✅ Next steps

---

## 🎬 What Happens During Deployment

```
[00:00] Installing Railway CLI...
[00:30] Logging in to Railway...
[01:00] Creating project...
[02:00] Provisioning PostgreSQL...
[03:00] Provisioning Redis...
[04:00] Setting environment variables...
[05:00] Deploying Backend API...
[10:00] Deploying Frontend...
[13:00] Deploying Worker...
[15:00] ✅ All services deployed!
```

**First deployment:** ~15 minutes (installing RDKit, PyTorch, etc.)
**Future deployments:** ~3-5 minutes (cached)

---

## 🎯 Alternative: Deploy via Railway Dashboard

If you prefer a visual interface:

1. **Go to**: https://railway.app/new
2. **Click**: "Deploy from GitHub repo"
3. **Select**: `MichaelCrowe11/crios-dr-crowe-coder`
4. **Branch**: `claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`
5. **Add Databases**:
   - Click "+ New" → "Database" → "PostgreSQL"
   - Click "+ New" → "Database" → "Redis"
6. **Set Variables**:
   - Go to each service → Variables
   - Add `ANTHROPIC_API_KEY=your_key`
   - Add `JWT_SECRET` (Railway can generate)
7. **Deploy**: Railway auto-deploys!

---

## 📊 After Deployment

### Check Status

```bash
railway status
```

You should see:
```
✓ crios-backend     Deployed    https://crios-backend-xxx.railway.app
✓ crios-frontend    Deployed    https://crios-frontend-xxx.railway.app
✓ crios-worker      Deployed
✓ postgres          Running
✓ redis             Running
```

### Test Your Deployment

```bash
# Test backend (replace with your URL from railway status)
curl https://crios-backend-xxx.railway.app/health

# Should return: {"status":"healthy"}

# Open frontend in browser
railway open crios-frontend
```

### View Logs

```bash
# Backend logs
railway logs --service crios-backend --follow

# Frontend logs
railway logs --service crios-frontend --follow
```

---

## 🐛 Common Issues

### "railway: command not found"

The script installs it automatically. If it fails:
```bash
npm i -g @railway/cli
```

### "Node.js not found"

Install Node.js from: https://nodejs.org
(Download the LTS version)

### "Permission denied: ./deploy-railway-cli.sh"

Make it executable:
```bash
chmod +x deploy-railway-cli.sh
```

### "Build failed"

Check logs:
```bash
railway logs --service crios-backend
```

Most common: Still building (wait a few minutes)

### "Frontend can't connect to backend"

```bash
# Get backend URL
railway status

# Set it for frontend
railway variables set NEXT_PUBLIC_API_URL=https://your-backend-url.railway.app --service crios-frontend

# Restart frontend
railway restart --service crios-frontend
```

---

## 💰 Cost

**Hobby Plan:** ~$10-15/month
- Includes $5 free credit (first month almost free!)
- 6 services (Backend, Frontend, Worker, PostgreSQL, Redis, Core API)
- Perfect for development/testing

**Pro Plan:** ~$25-35/month
- Better performance
- Custom domains
- Team features
- Automatic backups

---

## 🎓 What You're Deploying

### All Integrated Components ✅

- **Drug Discovery Engine**: RDKit-powered molecular analysis
- **ML Pipeline**: ChEMBL, PubMed, DrugBank integration
- **194 PhD Agent System**: AI orchestration for discoveries
- **Interactive Dashboard**: Real-time compound analysis
- **Background Processing**: Celery workers for heavy computation
- **Database**: PostgreSQL for compound storage
- **Caching**: Redis for performance
- **Authentication**: JWT-based secure access
- **Ethics & Safety**: Automated compliance screening

---

## 📚 Full Documentation

If you want to read more first:

| Document | Use Case |
|----------|----------|
| **`START_HERE_DEPLOY.md`** | **⭐ This file - simplest guide** |
| `RAILWAY_CLI_GUIDE.md` | Detailed CLI commands |
| `RAILWAY_DEPLOYMENT.md` | Complete 30-page guide |
| `DEPLOY_NOW.md` | Multiple deployment options |
| `QUICK_START_RAILWAY.md` | 5-minute overview |

---

## ✅ Deployment Checklist

- [ ] Railway account created
- [ ] Anthropic API key obtained
- [ ] Node.js 18+ installed (check: `node -v`)
- [ ] Repository cloned
- [ ] On correct branch (`claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`)
- [ ] Ran `./deploy-railway-cli.sh`
- [ ] Deployment completed successfully
- [ ] Services are healthy
- [ ] Frontend accessible

---

## 🆘 Need Help?

**Railway Support:**
- Discord: https://discord.gg/railway (very active!)
- Docs: https://docs.railway.app
- Status: https://status.railway.app

**CriOS Issues:**
- GitHub: https://github.com/MichaelCrowe11/crios-dr-crowe-coder/issues

---

## 🎉 Ready to Deploy!

### The Simplest Path:

```bash
cd ~/Desktop  # or wherever you want the code

git clone https://github.com/MichaelCrowe11/crios-dr-crowe-coder.git
cd crios-dr-crowe-coder
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH

./deploy-railway-cli.sh
```

**That's it! The script guides you through everything else.**

---

## ⏱️ Time Estimate

- Repository clone: 30 seconds
- Script prompts: 5 minutes
- **Deployment**: 10-15 minutes
- **Total**: ~20 minutes start to finish

---

**Your AI-powered drug discovery platform will be live in 20 minutes! 🚀💊🔬**

**Run the commands above on your computer and you're done!**

---

*All configuration is complete. All code is committed. Just run the script!* ✅
