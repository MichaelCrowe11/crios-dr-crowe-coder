# 🚀 CriOS - Ready to Deploy!

Everything is configured and ready for deployment. Follow these steps to deploy your application.

## ✅ What's Already Done

- ✅ Vercel CLI installed (v48.8.0)
- ✅ Vercel configuration files created
- ✅ Fly.io configuration files created
- ✅ Deployment documentation updated
- ✅ Helper scripts created
- ✅ All changes committed to git

## 🎯 Deployment Steps

### Step 1: Deploy to Vercel (Frontend)

#### Option A: Using Vercel CLI (Fastest)

```bash
# Login to Vercel
vercel login

# Deploy Frontend
cd platform/frontend
vercel --prod

# Deploy Algorithm Studio UI
cd ../../ui
vercel --prod
```

#### Option B: Using GitHub Integration (Recommended for Auto-Deploy)

1. **Go to Vercel**: https://vercel.com/new
2. **Import Repository**: `MichaelCrowe11/crios-dr-crowe-coder`
3. **Create First Project (Frontend)**:
   - Project Name: `crios-frontend`
   - Root Directory: `platform/frontend`
   - Framework: Next.js
   - Click "Deploy"

4. **Add Environment Variables** (after deployment):
   ```
   NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
   NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev
   NODE_ENV=production
   ```

5. **Create Second Project (UI)**:
   - Import same repository again
   - Project Name: `crios-ui`
   - Root Directory: `ui`
   - Framework: Next.js
   - Click "Deploy"

6. **Add Environment Variables**:
   ```
   NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
   NODE_ENV=production
   ```

#### Option C: Using Helper Script

```bash
./deploy-vercel.sh
```

---

### Step 2: Deploy Backend to Fly.io

```bash
# Install Fly.io CLI (if not already installed)
curl -L https://fly.io/install.sh | sh

# Login to Fly.io
flyctl auth login

# Create PostgreSQL database
flyctl postgres create --name crios-db --region iad --vm-size shared-cpu-1x --volume-size 10
# Save the DATABASE_URL from output

# Create Redis instance
flyctl redis create --name crios-redis --region iad --plan free
# Save the REDIS_URL from output

# Deploy Backend
cd platform/backend

# Create app (first time only)
flyctl apps create crios-backend --org personal

# Set secrets (replace with your actual values)
flyctl secrets set \
  DATABASE_URL="postgres://user:pass@host:5432/db" \
  REDIS_URL="redis://default:pass@host:6379" \
  ANTHROPIC_API_KEY="your-api-key-here"

# Deploy
flyctl deploy

# Check status
flyctl status
flyctl logs
```

---

## 🌐 Your Deployed URLs

After deployment, your application will be available at:

### Vercel (Frontend)
- **Frontend**: https://crios-frontend.vercel.app
- **Algorithm Studio**: https://crios-ui.vercel.app

### Fly.io (Backend)
- **API**: https://crios-backend.fly.dev
- **API Docs**: https://crios-backend.fly.dev/docs
- **Health Check**: https://crios-backend.fly.dev/health

---

## 🔑 Required Environment Variables

### For Vercel (Frontend)
```bash
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev
NODE_ENV=production
```

### For Vercel (UI)
```bash
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NODE_ENV=production
```

### For Fly.io (Backend)
```bash
DATABASE_URL=postgres://username:password@hostname:5432/database
REDIS_URL=redis://default:password@hostname:6379
ANTHROPIC_API_KEY=your-anthropic-api-key
PYTHONPATH=/app
ENV=production
```

---

## 🔍 Quick Deployment Commands

### Vercel (One-liner)
```bash
# Login once
vercel login

# Deploy both frontends
cd platform/frontend && vercel --prod && cd ../../ui && vercel --prod
```

### Fly.io (One-liner after setup)
```bash
cd platform/backend && flyctl deploy
```

---

## 🎨 Alternative: Deploy Everything via GitHub

### Vercel (Automatic)
1. Connect your GitHub repo to Vercel
2. Every push to `main` auto-deploys
3. Pull requests get preview deployments

### Fly.io (with GitHub Actions)

Create `.github/workflows/deploy-backend.yml`:

```yaml
name: Deploy Backend to Fly.io

on:
  push:
    branches: [main]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: superfly/flyctl-actions/setup-flyctl@master
      - run: flyctl deploy --remote-only -c platform/backend/fly.toml
        env:
          FLY_API_TOKEN: ${{ secrets.FLY_API_TOKEN }}
```

---

## 🛠 Troubleshooting

### Vercel Issues
```bash
# Check build locally
cd platform/frontend
npm install
npm run build

# View deployment logs
vercel logs crios-frontend
```

### Fly.io Issues
```bash
# View logs
flyctl logs -a crios-backend

# SSH into machine
flyctl ssh console -a crios-backend

# Check app status
flyctl status -a crios-backend

# Restart app
flyctl apps restart crios-backend
```

### Common Issues

1. **"No credentials found" on Vercel**
   - Run: `vercel login`

2. **Build fails on Vercel**
   - Check `package.json` dependencies
   - Ensure all env vars are set
   - Check build logs in Vercel dashboard

3. **Backend can't connect to DB**
   - Verify DATABASE_URL secret is set
   - Check PostgreSQL is running
   - Use `flyctl postgres db list`

4. **Frontend can't reach backend**
   - Verify NEXT_PUBLIC_API_URL is correct
   - Check CORS settings on backend
   - Ensure backend is deployed and running

---

## 📊 Deployment Checklist

- [ ] Install Vercel CLI: `npm i -g vercel`
- [ ] Install Fly.io CLI: `curl -L https://fly.io/install.sh | sh`
- [ ] Login to Vercel: `vercel login`
- [ ] Login to Fly.io: `flyctl auth login`
- [ ] Create PostgreSQL on Fly.io
- [ ] Create Redis on Fly.io
- [ ] Deploy backend to Fly.io
- [ ] Set backend secrets (DATABASE_URL, REDIS_URL, ANTHROPIC_API_KEY)
- [ ] Deploy frontend to Vercel
- [ ] Deploy UI to Vercel
- [ ] Set Vercel environment variables
- [ ] Test deployments
- [ ] Configure custom domains (optional)

---

## 💰 Cost Estimate

### Free Tier (Perfect for MVP)
- **Vercel**: Free (100GB bandwidth, 100 deployments/day)
- **Fly.io Backend**: Free (shared-cpu-1x, 256MB)
- **Fly.io PostgreSQL**: Free (1GB storage)
- **Fly.io Redis**: Free tier available
- **Total**: $0/month

### Production (Recommended)
- **Vercel Pro**: $20/month (1TB bandwidth)
- **Fly.io Backend**: ~$5-10/month (1GB RAM)
- **Fly.io PostgreSQL**: ~$5-10/month (10GB storage)
- **Fly.io Redis**: ~$5/month
- **Total**: ~$35-50/month

---

## 🎯 Quick Start (Fastest Way)

```bash
# 1. Deploy via Vercel GitHub integration (2 minutes)
# Visit: https://vercel.com/new
# Import: MichaelCrowe11/crios-dr-crowe-coder
# Create 2 projects (frontend + ui)

# 2. Deploy backend (5 minutes)
curl -L https://fly.io/install.sh | sh
flyctl auth login
flyctl postgres create --name crios-db --region iad
flyctl redis create --name crios-redis --region iad
cd platform/backend
flyctl apps create crios-backend
flyctl secrets set DATABASE_URL="..." REDIS_URL="..." ANTHROPIC_API_KEY="..."
flyctl deploy

# Done! 🎉
```

---

## 📚 Additional Resources

- Vercel Documentation: https://vercel.com/docs
- Fly.io Documentation: https://fly.io/docs
- Full Deployment Guide: See `DEPLOYMENT.md`
- Helper Script: Run `./deploy-vercel.sh`

---

## 🚨 Important Notes

1. **Backend must be deployed first** - Frontend needs backend URL
2. **Set environment variables** - Apps won't work without them
3. **ANTHROPIC_API_KEY** - Required for AI features
4. **Custom domains** - Configure after initial deployment
5. **Security vulnerabilities** - Run `npm audit fix` before production

---

## ✨ Next Steps After Deployment

1. ✅ Test all endpoints
2. ✅ Verify molecular visualization works
3. ✅ Test AI agent integrations
4. ✅ Run sample drug discovery workflows
5. ✅ Set up monitoring and alerting
6. ✅ Configure database backups
7. ✅ Set up custom domains
8. ✅ Enable CI/CD with GitHub Actions

---

**Ready to deploy? Start with Step 1 above!** 🚀

For issues or questions, see `DEPLOYMENT.md` for detailed troubleshooting.
